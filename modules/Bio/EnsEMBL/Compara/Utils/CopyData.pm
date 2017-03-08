=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2017] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=head1 NAME

Bio::EnsEMBL::Compara::Utils::CopyData

=head1 DESCRIPTION

This package exports method to copy_data between databases.
copy_data() is usually used to copy whole tables or large chunks of data,
without paying attention to the foreign-key constraints. It has two
specialized version copy_data_in_binary_mode() and copy_data_in_text_mode()
that are automatically used depending on the data-types in the table.
copy_data_with_foreign_keys_by_constraint() can copy individual rows with
their own depedencies. It will also "expand" the data, for instance by
copying homology_member too when asked to copy homology_member.

=head1 SYNOPSIS

  use Bio::EnsEMBL::Compara::Utils::CopyData (:row_copy);
  # The ", 1" at the end tells the function to "expand" the data, i.e. copy
  # extra rows to make the objects complete. Without it, it wouldn't copy
  # family_member
  copy_data_with_foreign_keys_by_constraint($source_dbc, $target_dbc, 'family', 'stable_id', 'ENSFM00730001521062', undef, 1);
  copy_data_with_foreign_keys_by_constraint($source_dbc, $target_dbc, 'gene_tree_root', 'stable_id', 'ENSGT00390000003602', undef, 1);

  # To insert a large number of rows in an optimal manner
  bulk_insert($self->compara_dba->dbc, 'homology_id_mapping', $self->param('homology_mapping'), ['mlss_id', 'prev_release_homology_id', 'curr_release_homology_id'], 'INSERT IGNORE');

=head1 AUTHORSHIP

Ensembl Team. Individual contributions can be found in the GIT log.

=cut

package Bio::EnsEMBL::Compara::Utils::CopyData;

use strict;
use warnings;

use base qw(Exporter);

our %EXPORT_TAGS;
our @EXPORT_OK;

@EXPORT_OK = qw(
    copy_data_with_foreign_keys_by_constraint
    clear_copy_data_cache
    copy_data
    copy_data_in_binary_mode
    copy_data_in_text_mode
    copy_table
    copy_table_in_binary_mode
    copy_table_in_text_mode
    bulk_insert
    single_insert
);
%EXPORT_TAGS = (
  'row_copy'    => [qw(copy_data_with_foreign_keys_by_constraint clear_copy_data_cache)],
  'table_copy'  => [qw(copy_data copy_table)],
  'insert'      => [qw(bulk_insert single_insert)],
  'all'         => [@EXPORT_OK]
);

use constant MAX_ROWS_FOR_MYSQLIMPORT => 1_000_000;

use Data::Dumper;
use File::Temp qw/tempfile/;

use Bio::EnsEMBL::Utils::Scalar qw(check_ref assert_ref);

my %foreign_key_cache = ();

my %data_expansions = (
    'ncbi_taxa_node' => [['taxon_id', 'ncbi_taxa_name', 'taxon_id'], ['parent_id', 'ncbi_taxa_node', 'taxon_id']],
    'gene_tree_root' => [['root_id', 'gene_tree_root_tag', 'root_id'], ['root_id', 'gene_tree_root_attr', 'root_id'], ['root_id', 'gene_tree_root', 'ref_root_id'], ['root_id', 'gene_tree_node', 'root_id'], ['root_id', 'homology', 'gene_tree_root_id'], ['root_id', 'CAFE_gene_family', 'gene_tree_root_id']],
    'CAFE_gene_family' => [['cafe_gene_family_id', 'CAFE_species_gene', 'cafe_gene_family_id']],
    'gene_tree_node' => [['node_id', 'gene_tree_node_tag', 'node_id'], ['node_id', 'gene_tree_node_attr', 'node_id']],
    'species_tree_node' => [['node_id', 'species_tree_node_tag', 'node_id'],['node_id', 'species_tree_node_attr', 'node_id'], ['root_id', 'species_tree_root', 'root_id'], ['parent_id', 'species_tree_node', 'node_id']],
    'species_tree_root' => [['root_id', 'species_tree_node', 'root_id']],
    'method_link_species_set' => [['method_link_species_set_id', 'method_link_species_set_tag', 'method_link_species_set_id'], ['method_link_species_set_id', 'method_link_species_set_attr', 'method_link_species_set_id']],
    'species_set_header' => [['species_set_id', 'species_set', 'species_set_id'], ['species_set_id', 'species_set_tag', 'species_set_id']],
    'family' => [['family_id', 'family_member', 'family_id']],
    'homology' => [['homology_id', 'homology_member', 'homology_id']],
    'gene_align' => [['gene_align_id', 'gene_align_member', 'gene_align_id']],
    'seq_member' => [['seq_member_id', 'other_member_sequence', 'seq_member_id']],
    'gene_member' => [['canonical_member_id', 'seq_member', 'seq_member_id']],
    'synteny_region' => [['synteny_region_id', 'dnafrag_region', 'synteny_region_id']],
    'genomic_align_block' => [['genomic_align_block_id', 'genomic_align', 'genomic_align_block_id'], ['genomic_align_block_id', 'conservation_score', 'genomic_align_block_id']],
    # FIXME: how to expand genomic_align_tree ?

);


=head2 copy_data_with_foreign_keys_by_constraint

  Arg[1]      : Bio::EnsEMBL::DBSQL::DBConnection $from_dbc
  Arg[2]      : Bio::EnsEMBL::DBSQL::DBConnection $to_dbc
  Arg[3]      : string $table_name
  Arg[4]      : (opt) string $where_field: the name of the column to use for the filtering
  Arg[5]      : (opt) string $where_value: the value of the column used for the filtering
  Arg[6]      : (opt) Bio::EnsEMBL::DBSQL::DBConnection $foreign_keys_dbc
  Arg[7]      : (bool) $expand_tables (default 0)

  Description : Copy the rows of this table that match the condition. Dependant rows
                (because of foreign keys) are automatically copied over. Foreign keys
                are discovered by asking the DBI layer. This is obviously only available
                for InnoDB, so if the schema is in MyISAM you'll have to provide
                $foreign_keys_dbc.
                If $expand_tables is set, the function will copy extra rows to make the
                objects "complete" (e.g. copy homology_member too if asked to copy homology)

=cut

sub copy_data_with_foreign_keys_by_constraint {
    my ($from_dbc, $to_dbc, $table_name, $where_field, $where_value, $foreign_keys_dbc, $expand_tables) = @_;
    assert_ref($from_dbc, 'Bio::EnsEMBL::DBSQL::DBConnection', 'from_dbc');
    assert_ref($to_dbc, 'Bio::EnsEMBL::DBSQL::DBConnection', 'to_dbc');
    die "A table name must be given" unless $table_name;
    my $fk_rules = _load_foreign_keys($foreign_keys_dbc, $to_dbc, $from_dbc);
    _memoized_insert($from_dbc, $to_dbc, $table_name, $where_field, $where_value, $fk_rules, $expand_tables);
}

# Load all the foreign keys and cache the result
sub _load_foreign_keys {

    foreach my $dbc (@_) {
        next unless check_ref($dbc, 'Bio::EnsEMBL::DBSQL::DBConnection');
        return $foreign_key_cache{$dbc->locator} if $foreign_key_cache{$dbc->locator};

        my $sth = $dbc->db_handle->foreign_key_info(undef, $dbc->dbname, undef, undef, undef, undef);
        my %fk = ();
        foreach my $x (@{ $sth->fetchall_arrayref() }) {
            # A.x REFERENCES B.y : push @{$fk{'A'}}, ['x', 'B', 'y'];
            push @{$fk{$x->[6]}}, [$x->[7], $x->[2], $x->[3]];
        }
        if (%fk) {
            $foreign_key_cache{$dbc->locator} = \%fk;
            return $foreign_key_cache{$dbc->locator};
        }
    }
    die "None of the DBConnections point to a database with foreign keys. Foreign keys are needed to populate linked tables\n";
}


my %cached_inserts = ();
sub _memoized_insert {
    my ($from_dbc, $to_dbc, $table, $where_field, $where_value, $fk_rules, $expand_tables) = @_;

    my $key = $table;
    $key .= "||$where_field=" . ($where_value // '__NA__') if $where_field;
    return if $cached_inserts{$from_dbc->locator}{$to_dbc->locator}{$key};
    $cached_inserts{$from_dbc->locator}{$to_dbc->locator}{$key} = 1;

    my $sql_select;
    my @execute_args;
    if ($where_field) {
        if (defined $where_value) {
            $sql_select = sprintf('SELECT * FROM %s WHERE %s = ?', $table, $where_field);
            push @execute_args, $where_value;
        } else {
            $sql_select = sprintf('SELECT * FROM %s WHERE %s IS NULL', $table, $where_field);
        }
    } else {
        $sql_select = 'SELECT * FROM '.$table;
    }
    #warn "<< $sql_select  using '@execute_args'\n";
    my $sth = $from_dbc->prepare($sql_select);
    $sth->execute(@execute_args);
    while (my $h = $sth->fetchrow_hashref()) {
        my %this_row = %$h;

        # First insert the requirements (to satisfy the foreign keys)
        _insert_related_rows($from_dbc, $to_dbc, \%this_row, $fk_rules->{$table}, $table, $where_field, $where_value, $fk_rules, $expand_tables);

        # Then the data
        my @cols = keys %this_row;
        my @qms  = map {'?'} @cols;
        my @vals = @this_row{@cols};
        my $insert_sql = sprintf('INSERT IGNORE INTO %s (%s) VALUES (%s)', $table, join(',', @cols), join(',', @qms));
        #warn ">> $insert_sql using '", join("','", map {$_//'<NULL>'} @vals), "'\n";
        my $rows = $to_dbc->do($insert_sql, undef, @vals);
        #warn "".($rows ? "true" : "false")." ".($rows == 0 ? "zero" : "non-zero")."\n";
        # no rows affected is translated into '0E0' which is true and ==0 at the same time
        if ($rows) {
            if ($rows == 0) {
                # The row was already there, but we can't assume it's been expanded too
            }
        } else {
            die "FAILED: ".$to_dbc->db_handle->errstr;
        }

        # And the expanded stuff
        _insert_related_rows($from_dbc, $to_dbc, \%this_row, $data_expansions{$table}, $table, $where_field, $where_value, $fk_rules, $expand_tables) if $expand_tables;
    }
}

sub _insert_related_rows {
    my ($from_dbc, $to_dbc, $this_row, $rules, $table, $where_field, $where_value, $fk_rules, $expand_tables) = @_;
    foreach my $x (@$rules) {
        #warn sprintf("%s(%s) needs %s(%s)\n", $table, @$x);
        if (not defined $this_row->{$x->[0]}) {
            next;
        } elsif (($table eq $x->[1]) and $where_field and ($where_field eq $x->[2]) and (defined $where_value) and ($where_value eq $this_row->{$x->[0]})) {
            # self-loop catcher: the code is about to store the same row again and again
            # we fall here when trying to insert a root gene_tree_node because its root_id links to itself
            next;
        }
        _memoized_insert($from_dbc, $to_dbc, $x->[1], $x->[2], $this_row->{$x->[0]}, $fk_rules, $expand_tables);
    }
}


=head2 clear_copy_data_cache

  Description : Clears the cache to free up some memory

=cut

sub clear_copy_data_cache {
    %cached_inserts = ();
}


=head2 _has_binary_column

  Example     : _has_binary_column($dbc, 'genomic_align_block');
  Description : Tells whether the table has a binary column
  Returntype  : Boolean
  Exceptions  : none
  Caller      : general
  Status      : Stable

=cut

sub _has_binary_column {
    my ($dbc, $table_name) = @_;

    assert_ref($dbc, 'Bio::EnsEMBL::DBSQL::DBConnection', 'dbc');
    return $dbc->{"_has_binary_column__${table_name}"} if exists $dbc->{"_has_binary_column__${table_name}"};

    my $sth = $dbc->db_handle->column_info($dbc->dbname, undef, $table_name, '%');
    $sth->execute;
    my $all_rows = $sth->fetchall_arrayref;
    my $binary_mode = 0;
    foreach my $this_col (@$all_rows) {
        if (($this_col->[5] =~ /BINARY$/) or ($this_col->[5] =~ /BLOB$/) or ($this_col->[5] eq "BIT")) {
            $binary_mode = 1;
            last;
        }
    }
    $dbc->{"_has_binary_column__${table_name}"} = $binary_mode;
    return $binary_mode;
}


=head2 copy_data

  Arg[1]      : Bio::EnsEMBL::DBSQL::DBConnection $from_dbc
  Arg[2]      : Bio::EnsEMBL::DBSQL::DBConnection $to_dbc
  Arg[3]      : string $table_name
  Arg[4]      : (opt) string $query
  Arg[5]      : (opt) string $index_name
  Arg[6]      : (opt) integer $min_id
  Arg[7]      : (opt) integer $max_id
  Arg[8]      : (opt) integer $step
  Arg[9]      : (opt) boolean $disable_keys (default: true)
  Arg[10]     : (opt) boolean $reenable_keys (default: true)
  Arg[11]     : (opt) boolean $holes_possible (default: false)
  Arg[12]     : (opt) boolean $replace (default: false) [only used when the underlying data is text-only]

  Description : Copy data in this table. The main optional arguments are:
                 - ($index_name,$min_id,$max_id) to restrict to a range
                 - $query for a general way of altering the data (e.g. fixing it)

=cut

sub copy_data {
    my ($from_dbc, $to_dbc, $table_name, $query, $index_name, $min_id, $max_id, $step, $disable_keys, $reenable_keys, $holes_possible, $replace) = @_;

    assert_ref($from_dbc, 'Bio::EnsEMBL::DBSQL::DBConnection', 'from_dbc');
    assert_ref($to_dbc, 'Bio::EnsEMBL::DBSQL::DBConnection', 'to_dbc');

    print "Copying data in table $table_name\n";

    die "Keys must be disabled in order to be reenabled\n" if $reenable_keys and not $disable_keys;

    #When merging the patches, there are so few items to be added, we have no need to disable the keys
    if ($disable_keys // 1) {
        #speed up writing of data by disabling keys, write the data, then enable
        #but takes far too long to ENABLE again
        $to_dbc->do("ALTER TABLE `$table_name` DISABLE KEYS");
    }
    if (_has_binary_column($from_dbc, $table_name)) {
        copy_data_in_binary_mode($from_dbc, $to_dbc, $table_name, $query, $index_name, $min_id, $max_id, $step);
    } else {
        copy_data_in_text_mode($from_dbc, $to_dbc, $table_name, $query, $index_name, $min_id, $max_id, $step, $holes_possible, $replace);
    }
    if ($reenable_keys // 1) {
        $to_dbc->do("ALTER TABLE `$table_name` ENABLE KEYS");
    }
}


=head2 _escape

  Description : Helper function that escapes some special characters.

=cut

sub _escape {
    my $s = shift;
    return '\N' unless defined $s;
    $s =~ s/\n/\\\n/g;
    $s =~ s/\t/\\\t/g;
    return $s;
}


=head2 copy_data_in_text_mode

  Description : A specialized version of copy_data() for tables that don't have
                any binary data and can be loaded with mysqlimport.

=cut

sub copy_data_in_text_mode {
    my ($from_dbc, $to_dbc, $table_name, $query, $index_name, $min_id, $max_id, $step, $holes_possible, $replace) = @_;

    assert_ref($from_dbc, 'Bio::EnsEMBL::DBSQL::DBConnection', 'from_dbc');
    assert_ref($to_dbc, 'Bio::EnsEMBL::DBSQL::DBConnection', 'to_dbc');

    my $user = $to_dbc->username;
    my $pass = $to_dbc->password;
    my $host = $to_dbc->host;
    my $port = $to_dbc->port;
    my $dbname = $to_dbc->dbname;

    #Default step size.
    $step ||= MAX_ROWS_FOR_MYSQLIMPORT;

    my ($use_limit, $start);
    if (defined $index_name && defined $min_id && defined $max_id) {
        # We'll use BETWEEN
        $use_limit = 0;
        $start = $min_id;
    } else {
        # We'll use LIMIT 
        $use_limit = 1;
        $start = 0;
    }
    # $use_limit also tells whether $start and $end are counters or values comparable to $index_name

    my $total_rows = 0;
    while (1) {
        #my $start_time = time();
        my $end = $start + $step - 1;
        my $sth;
        my $sth_attribs = { 'mysql_use_result' => 1 };

        #print "start $start end $end\n";
        #print "query $query\n";
        if ($use_limit) {
            $sth = $from_dbc->prepare( $query." LIMIT $start, $step", $sth_attribs );
        } else {
            $sth = $from_dbc->prepare( $query." AND $index_name BETWEEN $start AND $end", $sth_attribs );
        }
        $start += $step;
        $sth->execute();
        my $first_row = $sth->fetchrow_arrayref;

        ## EXIT CONDITION
        if (!$first_row) {
            # We're told there could be holes in the data, and $end hasn't yet reached $max_id
            next if ($holes_possible and !$use_limit and ($end < $max_id));
            # Otherwise it is the end
            return $total_rows;
        }

        #my $time = time();
        my ($fh, $filename) = tempfile("${table_name}.XXXXXX", TMPDIR => 1);
        print $fh join("\t", map {_escape($_)} @$first_row), "\n";
        my $nrows = 1;
        while(my $this_row = $sth->fetchrow_arrayref) {
            print $fh join("\t", map {_escape($_)} @$this_row), "\n";
            $nrows++;
        }
        close($fh);
        $sth->finish;
        #print "start $start end $end $max_id rows $nrows\n";
        #print "FILE $filename\n";
        #print "time " . ($start-$min_id) . " " . (time - $start_time) . "\n";

        # Heuristics: it's going to take some time to insert these rows, so
        # better to disconnect to save resources on the source server
        $from_dbc->disconnect_if_idle if $nrows >= 100_000;

        system('mysqlimport', "-h$host", "-P$port", "-u$user", $pass ? ("-p$pass") : (), '--local', '--lock-tables', $replace ? '--replace' : '--ignore', $dbname, $filename);

        unlink($filename);
        $total_rows += $nrows;
    }
    return $total_rows;
}


=head2 copy_data_in_binary_mode

  Description : A specialized version of copy_data() for tables that have binary
                data, using mysqldump.

=cut

sub copy_data_in_binary_mode {
    my ($from_dbc, $to_dbc, $table_name, $query, $index_name, $min_id, $max_id, $step) = @_;

    assert_ref($from_dbc, 'Bio::EnsEMBL::DBSQL::DBConnection', 'from_dbc');
    assert_ref($to_dbc, 'Bio::EnsEMBL::DBSQL::DBConnection', 'to_dbc');

    my $from_user = $from_dbc->username;
    my $from_pass = $from_dbc->password;
    my $from_host = $from_dbc->host;
    my $from_port = $from_dbc->port;
    my $from_dbname = $from_dbc->dbname;

    my $to_user = $to_dbc->username;
    my $to_pass = $to_dbc->password;
    my $to_host = $to_dbc->host;
    my $to_port = $to_dbc->port;
    my $to_dbname = $to_dbc->dbname;

    my $use_limit = 0;
    my $start = $min_id;

    #all the data in the table needs to be copied and does not need fixing
    if (!defined $query) {
        #my $start_time  = time();

        system("mysqldump -h$from_host -P$from_port -u$from_user ".($from_pass ? "-p$from_pass" : '')." --insert-ignore -t $from_dbname $table_name ".
            "| mysql   -h$to_host   -P$to_port   -u$to_user   ".($to_pass ? "-p$to_pass" : '')." $to_dbname");

        #print "time " . ($start-$min_id) . " " . (time - $start_time) . "\n";

        return;
    }

    print " ** WARNING ** Copying table $table_name in binary mode, this requires write access.\n";
    print " ** WARNING ** The original table will be temporarily renamed as original_$table_name.\n";
    print " ** WARNING ** An auxiliary table named temp_$table_name will also be created.\n";
    print " ** WARNING ** You may have to undo this manually if the process crashes.\n\n";

    #If not using BETWEEN, revert back to LIMIT
    if (!defined $index_name && !defined $min_id && !defined $max_id) {
        $use_limit = 1;
        $start = 0;
    }
    #my $start = 0;
    if (!defined $step) {
        $step = 1000000;
    }
    while (1) {
        #my $start_time  = time();
        my $end = $start + $step - 1;
        #print "start $start end $end\n";

        ## Copy data into a aux. table
        my $sth;
        if (!$use_limit) {
            $sth = $from_dbc->prepare("CREATE TABLE temp_$table_name $query AND $index_name BETWEEN $start AND $end");
        } else {
            $sth = $from_dbc->prepare("CREATE TABLE temp_$table_name $query LIMIT $start, $step");
        }
        $sth->execute();

        $start += $step;
        my $count = $from_dbc->db_handle->selectrow_array("SELECT count(*) FROM temp_$table_name");

        ## EXIT CONDITION
        if (!$count) {
            $from_dbc->db_handle->do("DROP TABLE temp_$table_name");
            return;
        }

        ## Change table names (mysqldump will keep the table name, hence we need to do this)
        $from_dbc->db_handle->do("ALTER TABLE $table_name RENAME original_$table_name");
        $from_dbc->db_handle->do("ALTER TABLE temp_$table_name RENAME $table_name");

        ## mysqldump data

        system("mysqldump -h$from_host -P$from_port -u$from_user ".($from_pass ? "-p$from_pass" : '')." --insert-ignore -t $from_dbname $table_name ".
            "| mysql   -h$to_host   -P$to_port   -u$to_user   ".($to_pass ? "-p$to_pass" : '')." $to_dbname");

        #print "time " . ($start-$min_id) . " " . (time - $start_time) . "\n";

        ## Undo table names change
        $from_dbc->db_handle->do("DROP TABLE $table_name");
        $from_dbc->db_handle->do("ALTER TABLE original_$table_name RENAME $table_name");

        #print "total time " . ($start-$min_id) . " " . (time - $start_time) . "\n";
    }
}


=head2 copy_table

  Arg[1]      : Bio::EnsEMBL::DBSQL::DBConnection $from_dbc
  Arg[2]      : Bio::EnsEMBL::DBSQL::DBConnection $to_dbc
  Arg[3]      : string $table_name
  Arg[4]      : (opt) string $where_filter
  Arg[5]      : (opt) boolean $replace (default: false)

  Description : Copy the table (either all of it or a subset).
                The main optional argument is $where_filter, which allows to select a portion of
                the table. Note: the filter must be valid on the table alone, and does not support
                JOINs. If you need the latter, use copy_data()

=cut

sub copy_table {
    my ($from_dbc, $to_dbc, $table_name, $where_filter, $replace, $skip_disable_keys) = @_;

    assert_ref($from_dbc, 'Bio::EnsEMBL::DBSQL::DBConnection', 'from_dbc');
    assert_ref($to_dbc, 'Bio::EnsEMBL::DBSQL::DBConnection', 'to_dbc');

    print "Copying data in table $table_name\n";

    if (_has_binary_column($from_dbc, $table_name)) {
        return copy_table_in_binary_mode($from_dbc, $to_dbc, $table_name, $where_filter, $replace, $skip_disable_keys);
    } else {
        return copy_table_in_text_mode($from_dbc, $to_dbc, $table_name, $where_filter, $replace);
    }
}


=head2 copy_table_in_text_mode

  Description : A specialized version of copy_table() for tables that don't have
                any binary data and can be loaded with mysqlimport.

=cut

sub copy_table_in_text_mode {
    my ($from_dbc, $to_dbc, $table_name, $where_filter, $replace) = @_;

    my $query = 'SELECT * FROM '.$table_name.($where_filter ? ' WHERE '.$where_filter : '');
    return copy_data_in_text_mode($from_dbc, $to_dbc, $table_name, $query, undef, undef, undef, undef, undef, $replace);
}


=head2 copy_table_in_binary_mode

  Description : A specialized version of copy_table() for tables that have binary
                data, using mysqldump.

=cut

sub copy_table_in_binary_mode {
    my ($from_dbc, $to_dbc, $table_name, $where_filter, $replace, $skip_disable_keys) = @_;

    assert_ref($from_dbc, 'Bio::EnsEMBL::DBSQL::DBConnection', 'from_dbc');
    assert_ref($to_dbc, 'Bio::EnsEMBL::DBSQL::DBConnection', 'to_dbc');

    my $from_user = $from_dbc->username;
    my $from_pass = $from_dbc->password;
    my $from_host = $from_dbc->host;
    my $from_port = $from_dbc->port;
    my $from_dbname = $from_dbc->dbname;

    my $to_user = $to_dbc->username;
    my $to_pass = $to_dbc->password;
    my $to_host = $to_dbc->host;
    my $to_port = $to_dbc->port;
    my $to_dbname = $to_dbc->dbname;

    #my $start_time  = time();
    my $insert_mode = $replace ? '--replace' : '--insert-ignore';

    system("mysqldump -h$from_host -P$from_port -u$from_user ".($from_pass ? "-p$from_pass" : '')." $insert_mode -t $from_dbname $table_name ".
        ($where_filter ? "-w '$where_filter'" : "")." ".
        ($skip_disable_keys ? "--skip-disable-keys" : "")." ".
        "| mysql   -h$to_host   -P$to_port   -u$to_user   ".($to_pass ? "-p$to_pass" : '')." $to_dbname");

    #print "time " . ($start-$min_id) . " " . (time - $start_time) . "\n";
}


=head2 single_insert

  Arg[1]      : Bio::EnsEMBL::DBSQL::DBConnection $dest_dbc
  Arg[2]      : string $table_name
  Arg[3]      : arrayref of values$data
  Arg[4]      : (opt) arrayref of strings $col_names (defaults to the column-order at the database level)
  Arg[5]      : (opt) string $insertion_mode (default: 'INSERT')

  Description : Simple method to execute an INSERT statement without having to write it. The values
                in $data must be in the same order as in $col_names (if provided) or the columns in
                the table itself.  The method returns the number of rows inserted (0 or 1).
  Returntype  : integer
  Exceptions  : none
  Caller      : general
  Status      : Stable

=cut

sub single_insert {
    my ($dest_dbc, $table_name, $data, $col_names, $insertion_mode) = @_;

    my $n_values = scalar(@$data);
    my $insert_sql = ($insertion_mode || 'INSERT') . ' INTO ' . $table_name;
    $insert_sql .= ' (' . join(',', @$col_names) . ')' if $col_names;
    $insert_sql .= ' VALUES (' . ('?,'x($n_values-1)) . '?)';
    return $dest_dbc->do($insert_sql, undef, @$data) or die "Could not execute the insert because of ".$dest_dbc->db_handle->errstr;
}


=head2 bulk_insert

  Arg[1]      : Bio::EnsEMBL::DBSQL::DBConnection $dest_dbc
  Arg[2]      : string $table_name
  Arg[3]      : arrayref of arrayrefs $data
  Arg[4]      : (opt) arrayref of strings $col_names (defaults to the column-order at the database level)
  Arg[5]      : (opt) string $insertion_mode (default: 'INSERT')

  Description : Execute extended INSERT statements (or whatever flavour selected in $insertion_mode)
                on $dest_dbc to push the data kept in the arrayref $data.  Each arrayref in $data
                corresponds to a row, in the same order as in $col_names (if provided) or the columns
                in the table itself.  The method returns the total number of rows inserted
  Returntype  : integer
  Exceptions  : none
  Caller      : general
  Status      : Stable

=cut

sub bulk_insert {
    my ($dest_dbc, $table_name, $data, $col_names, $insertion_mode) = @_;

    my $insert_n   = 0;
    while (@$data) {
        my $insert_sql = ($insertion_mode || 'INSERT') . ' INTO ' . $table_name;
        $insert_sql .= ' (' . join(',', @$col_names) . ')' if $col_names;
        $insert_sql .= ' VALUES ';
        my $first = 1;
        while (@$data and (length($insert_sql) < 1_000_000)) {
            $insert_sql .= ($first ? '' : ', ') . '(' . join(',', map {defined $_ ? '"'.$_.'"' : 'NULL'} @{shift @$data}) . ')';
            $first = 0;
        }
        $insert_n += $dest_dbc->do($insert_sql) or die "Could not execute the insert because of ".$dest_dbc->db_handle->errstr;
    }
    return $insert_n;
}


1;
 
