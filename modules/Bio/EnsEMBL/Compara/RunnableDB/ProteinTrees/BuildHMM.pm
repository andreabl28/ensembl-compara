=head1 LICENSE

  Copyright (c) 1999-2012 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

   http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=head1 NAME

Bio::EnsEMBL::Compara::RunnableDB::ProteinTrees::BuildHMM

=head1 DESCRIPTION

This Analysis/RunnableDB is designed to take a ProteinTree as input
This must already have a multiple alignment run on it. It uses that alignment
as input create a HMMER HMM profile

input_id/parameters format eg: "{'protein_tree_id'=>1234}"
    protein_tree_id : use 'id' to fetch a cluster from the ProteinTree

=head1 SYNOPSIS

my $db           = Bio::EnsEMBL::Compara::DBAdaptor->new($locator);
my $buildhmm = Bio::EnsEMBL::Compara::RunnableDB::ProteinTrees::BuildHMM->new
  (
   -db         => $db,
   -input_id   => $input_id,
   -analysis   => $analysis
  );
$buildhmm->fetch_input(); #reads from DB
$buildhmm->run();
$buildhmm->output();
$buildhmm->write_output(); #writes to DB

=head1 AUTHORSHIP

Ensembl Team. Individual contributions can be found in the CVS log.

=head1 MAINTAINER

$Author$

=head VERSION

$Revision$

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with an underscore (_)

=cut

package Bio::EnsEMBL::Compara::RunnableDB::ProteinTrees::BuildHMM;

use strict;

use IO::File;
use File::Basename;
use Time::HiRes qw(time gettimeofday tv_interval);

use Bio::AlignIO;
use Bio::SimpleAlign;

use Bio::EnsEMBL::Compara::Member;
use Bio::EnsEMBL::Compara::Graph::NewickParser;

use base ('Bio::EnsEMBL::Compara::RunnableDB::BaseRunnable');

sub param_defaults {
    return {
            'cdna'                  => 0,
    };
}


sub fetch_input {
    my $self = shift @_;

    $self->check_if_exit_cleanly;

    my $protein_tree_id     = $self->param('protein_tree_id') or die "'protein_tree_id' is an obligatory parameter";
    my $protein_tree        = $self->compara_dba->get_ProteinTreeAdaptor->fetch_node_by_node_id( $protein_tree_id )
                                        or die "Could not fetch protein_tree with protein_tree_id='$protein_tree_id'";
    $self->param('protein_tree', $protein_tree);

    my $hmm_type = $self->param('cdna') ? 'dna' : 'aa';

    if ($self->param('notaxon')) {
        $hmm_type .= "_notaxon" . "_" . $self->param('notaxon');
    }
    if ($self->param('taxon_ids')) {
        $hmm_type .= "_" . join(':', @{$self->param('taxon_ids')});
    }
    $self->param('hmm_type', $hmm_type);

    my $node_id = $self->param('protein_tree')->node_id;
    my $table_name = 'protein_tree_hmmprofile';
    my $query = "SELECT hmmprofile FROM $table_name WHERE type=\"$hmm_type\" AND node_id=$node_id";
    print STDERR "$query\n" if ($self->debug);
    my $sth = $self->compara_dba->dbc->prepare($query);
    $sth->execute;
    my $result = $sth->fetch;
    if (defined($result)) { # Has been done already
        $self->param('done', 1);
        return;
    }

  my @to_delete;

  if ($self->param('notaxon')) {
    foreach my $leaf (@{$protein_tree->get_all_leaves}) {
      next unless ($leaf->taxon_id eq $self->param('notaxon'));
      push @to_delete, $leaf;
    }
    $protein_tree = $protein_tree->remove_nodes(\@to_delete);
  }

  if ($self->param('taxon_ids')) {
    my $taxon_ids_to_keep;
    foreach my $taxon_id (@{$self->param('taxon_ids')}) {
      $taxon_ids_to_keep->{$taxon_id} = 1;
    }
    foreach my $leaf (@{$protein_tree->get_all_leaves}) {
      next if (defined($taxon_ids_to_keep->{$leaf->taxon_id}));
      push @to_delete, $leaf;
    }
    $protein_tree = $protein_tree->remove_nodes(\@to_delete);
  }

  if (!defined($protein_tree)) {
    $self->param('done', 1);
  }

  if (2 > (scalar @{$protein_tree->get_all_leaves})) {
    $self->param('done', 1);
  }
}


=head2 run

    Title   :   run
    Usage   :   $self->run
    Function:   runs hmmbuild
    Returns :   none
    Args    :   none

=cut


sub run {
    my $self = shift @_;

    $self->check_if_exit_cleanly;

    unless($self->param('done')) {
        $self->run_buildhmm;
    }
}


=head2 write_output

    Title   :   write_output
    Usage   :   $self->write_output
    Function:   stores hmmprofile
    Returns :   none
    Args    :   none

=cut


sub write_output {
    my $self = shift @_;

    $self->check_if_exit_cleanly;

    unless($self->param('done')) {
        $self->store_hmmprofile;
    }
}


sub DESTROY {
  my $self = shift;

  if($self->param('protein_tree')) {
    printf("BuildHMM::DESTROY  releasing tree\n") if($self->debug);
    $self->param('protein_tree')->release_tree;
    $self->param('protein_tree', undef);
  }

  $self->SUPER::DESTROY if $self->can("SUPER::DESTROY");
}


##########################################
#
# internal methods
#
##########################################


sub run_buildhmm {
  my $self = shift;

  my $starttime = time()*1000;

  my $stk_file = $self->dumpTreeMultipleAlignmentToWorkdir ( $self->param('protein_tree') ) or return;
  return if($self->param('done'));

  my $hmm_file = $self->param('hmm_file', $stk_file . '_hmmbuild.hmm');

  my $buildhmm_exe = $self->param('buildhmm_exe')
        or die "'buildhmm_exe' is an obligatory parameter";

  die "Cannot execute '$buildhmm_exe'" unless(-x $buildhmm_exe);

  ## as in treefam
  # $hmmbuild --amino -g -F $file.hmm $file >/dev/null

  my $cmd = $buildhmm_exe;
  $cmd .= ($self->param('cdna') ? ' --dna ' : ' --amino ');

  $cmd .= $hmm_file;
  $cmd .= " ". $stk_file;
  $cmd .= " 2>&1 > /dev/null" unless($self->debug);

  $self->compara_dba->dbc->disconnect_when_inactive(1);
  print("$cmd\n") if($self->debug);
  my $worker_temp_directory = $self->worker_temp_directory;
  $cmd = "cd $worker_temp_directory ; $cmd";
  if(system($cmd)) {
    my $system_error = $!;
    die "Could not run [$cmd] : $system_error";
  }

  $self->compara_dba->dbc->disconnect_when_inactive(0);
  my $runtime = time()*1000-$starttime;

  $self->param('protein_tree')->tree->store_tag('BuildHMM_runtime_msec', $runtime);
}


########################################################
#
# ProteinTree input/output section
#
########################################################

sub dumpTreeMultipleAlignmentToWorkdir {
  my $self = shift;
  my $protein_tree = shift;
  
  my $leafcount = scalar(@{$protein_tree->get_all_leaves});

  my $file_root = $self->worker_temp_directory. $protein_tree->node_id;
  $file_root =~ s/\/\//\//g;  # converts any // in path to /

  my $aln_file = $file_root . '.aln';
#  return $aln_file if(-e $aln_file);
  if($self->debug) {
    printf("dumpTreeMultipleAlignmentToWorkdir : %d members\n", $leafcount);
    print("aln_file = '$aln_file'\n");
  }

  open(OUTSEQ, ">$aln_file") or die "Could not open '$aln_file' for writing : $!";

  my $sa = $protein_tree->get_SimpleAlign (
     -id_type => 'MEMBER',
     -cdna => $self->param('cdna'),
     -stop2x => 1
  );
  $sa->set_displayname_flat(1);

  # Pairwise alns can sometimes be empty
  if (0 == scalar($sa->each_seq)) {
    return $self->param('done', 1);
  }

  my $alignIO = Bio::AlignIO->newFh
    (
     -fh => \*OUTSEQ,
     -format => "fasta"
    );
  print $alignIO $sa;

  close OUTSEQ;

  unless(-e $aln_file and -s $aln_file) {
    die "There are no alignments in '$aln_file', cannot continue";
  }

  my $stk_file = $file_root . '.stk';

  my $sreformat_exe = $self->param('sreformat_exe')
        or die "'sreformat_exe' is an obligatory parameter";

  die "Cannot execute '$sreformat_exe'" unless(-x $sreformat_exe);

  my $cmd = "$sreformat_exe stockholm $aln_file > $stk_file";
  if(system($cmd)) {
    my $system_error = $!;
    die "Could not run [$cmd] : $system_error";
  }
  unless(-e $stk_file and -s $stk_file) {
    die "'$cmd' did not produce any data in '$stk_file'";
  }

  return $stk_file;
}


sub store_hmmprofile {
  my $self = shift;
  my $hmm_file =  $self->param('hmm_file');
  my $protein_tree = $self->param('protein_tree');
  
  #parse hmmer file
  print("load from file $hmm_file\n") if($self->debug);
  open(FH, $hmm_file) or die "Could not open '$hmm_file' for reading : $!";
  my $hmm_text = join('', <FH>);
  close(FH);

  my $table_name = 'protein_tree_hmmprofile';
  my $sth = $self->compara_dba->dbc->prepare("INSERT INTO $table_name VALUES (?,?,?)");
  $sth->execute($protein_tree->node_id, $self->param('hmm_type'), $hmm_text);
}

1;
