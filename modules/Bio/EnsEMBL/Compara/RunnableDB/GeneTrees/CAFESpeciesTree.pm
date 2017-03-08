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

=cut

=head1 NAME

Bio::EnsEMBL::Compara::RunnableDB::GeneTrees::CAFESpeciesTree

=head1 DESCRIPTION

This RunnableDB builds a CAFE-compliant species tree (binary & ultrametric with time units).

=head1 INHERITANCE TREE

Bio::EnsEMBL::Compara::RunnableDB::BaseRunnable

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with an underscore (_)

=cut

package Bio::EnsEMBL::Compara::RunnableDB::GeneTrees::CAFESpeciesTree;

use strict;
use warnings;
use Data::Dumper;
use Scalar::Util qw(looks_like_number);
use Bio::EnsEMBL::Compara::Graph::NewickParser;

use base ('Bio::EnsEMBL::Compara::RunnableDB::BaseRunnable');

sub param_defaults {
    return {
            'tree_fmt'   => '%{n}%{":"d}',
            'label'      => 'full_species_tree',
            'new_label'  => 'cafe',
           };
}


=head2 fetch_input

    Title     : fetch_input
    Usage     : $self->fetch_input
    Function  : Fetches input data from database
    Returns   : none
    Args      : none

=cut

sub fetch_input {
    my ($self) = @_;

    $self->param_required('mlss_id');

    my $speciesTree_Adaptor = $self->compara_dba->get_SpeciesTreeAdaptor();
    $self->param('speciesTree_Adaptor', $speciesTree_Adaptor);

    my $genomeDB_Adaptor = $self->compara_dba->get_GenomeDBAdaptor();
    $self->param('genomeDB_Adaptor', $genomeDB_Adaptor);

    my $NCBItaxon_Adaptor = $self->compara_dba->get_NCBITaxon(); # Adaptor??
    $self->param('NCBItaxon_Adaptor', $NCBItaxon_Adaptor);

    my $CAFETree_Adaptor = $self->compara_dba->get_CAFEGeneFamilyAdaptor();
    $self->param('cafeTree_Adaptor', $CAFETree_Adaptor);

    my $full_species_tree = $self->compara_dba->get_SpeciesTreeAdaptor->fetch_by_method_link_species_set_id_label($self->param('mlss_id'), $self->param('label'));
    $self->param('full_species_tree', $full_species_tree); ## This is the full tree, not the string

    my $cafe_species = $self->param('cafe_species') || [];
    if (not ref($cafe_species)) {
        my $cafe_species_str = $self->param('cafe_species');
        $cafe_species_str =~ s/["'\[\] ]//g;
        $cafe_species = [split(',', $cafe_species_str)];
    }
    if (scalar(@{$cafe_species}) == 0) {  # No species for the tree. Make a full tree
        print STDERR "No species provided for the CAFE tree. I will take them all\n" if ($self->debug());
        $self->param('cafe_species', undef);
        $self->param('n_missing_species_in_tree', 0);
    } else {
        my %gdb_ids = map {$_->dbID => 1} map {$genomeDB_Adaptor->fetch_by_name_assembly($_) || die "Could not find a GenomeDB named '$_'"} @$cafe_species;
        $self->param('cafe_species', \%gdb_ids);
        $self->param('n_missing_species_in_tree', scalar(@{$genomeDB_Adaptor->fetch_all()})-scalar(@{$cafe_species}));
    }

    return;
}

sub run {
    my ($self) = @_;
    my $species_tree = $self->param('full_species_tree');
    my $species_tree_root = $species_tree->root;
    $species_tree_root->print_tree(0.2);
    my $species = $self->param('cafe_species');
    my $fmt = $self->param('tree_fmt');
    my $mlss_id = $self->param('mlss_id');
    print STDERR Dumper $species if ($self->debug());
#    my $eval_species_tree = Bio::EnsEMBL::Compara::Graph::NewickParser::parse_newick_into_tree($species_tree_string);

    $self->include_distance_to_parent($species_tree_root);
    $self->fix_ensembl_timetree_mya($species_tree_root);
    $self->ensembl_timetree_mya_to_distance_to_parent($species_tree_root);
    $self->ultrametrize($species_tree_root);
    my $binTree = $self->binarize($species_tree_root);
    $self->fix_zeros($binTree);

    my $cafe_tree_root;
    if (defined $species) {
        $cafe_tree_root = $self->prune_tree($binTree, $species);
    } else {
        $cafe_tree_root = $binTree;
    }
    $cafe_tree_root->distance_to_parent(0); # NULL would be more accurate
    $self->check_tree($cafe_tree_root);
    $cafe_tree_root->build_leftright_indexing();

    ## The modified tree is put back in the species tree object
    $species_tree->root($cafe_tree_root);

    # Store the tree (At this point, it is a species tree not a CAFE tree)
    my $speciesTree_Adaptor = $self->param('speciesTree_Adaptor');

    my $cafe_tree_str = $cafe_tree_root->newick_format('ryo', $fmt);
    print STDERR "Tree to store:\n$cafe_tree_str\n" if ($self->debug);

    $species_tree->label($self->param_required('new_label'));
    $speciesTree_Adaptor->store($species_tree);

}

sub write_output {
    my ($self) = @_;
    $self->dataflow_output_id( {
        'species_tree_root_id' => $self->param('full_species_tree')->root_id,
        'n_missing_species_in_tree' => $self->param('n_missing_species_in_tree'),
    }, 2);
}


#############################
## Internal methods #########
#############################



sub include_distance_to_parent {
    my ($self, $tree) = @_;
    my $NCBItaxon_Adaptor = $self->param('NCBItaxon_Adaptor');
    my $nodes = $tree->get_all_nodes();
    for my $node (@$nodes) {
        unless ($node->is_leaf) {
            my $taxon_id = $node->taxon_id();
            my $ncbiTaxon = $NCBItaxon_Adaptor->fetch_node_by_taxon_id($taxon_id);
            my $mya = $ncbiTaxon->get_value_for_tag('ensembl timetree mya');
            if (!$mya && $self->param('use_timetree_times')) {
                die "No 'ensembl timetree mya' tag for taxon_id=$taxon_id\n";
            }
            for my $child (@{$node->children}) {
                $child->distance_to_parent(int($mya));
            }
        }
    }
}

sub fix_ensembl_timetree_mya {
    my ($self, $tree) = @_;
    my $leaves = $tree->get_all_leaves();
    for my $leaf (@$leaves) {
        fix_path($leaf);
    }
}

sub fix_path {
    my ($node) = @_;
    for (;;) {
        if ($node->has_parent()) {
            if ($node->parent->distance_to_parent() == 0) {
                $node = $node->parent;
                next;
            }
            if ($node->parent()->distance_to_parent() < $node->distance_to_parent()) {
                # The if is because the root doesn't have proper mya set
                if ($node->parent->has_parent) {
                    die "'ensembl timetree mya' tags are not monotonous. Check ".$node->taxon_id." and ".$node->parent->taxon_id."\n";
                }
            }
        } else {
            return
        }
        $node = $node->parent();
    }
}

sub ensembl_timetree_mya_to_distance_to_parent {
    my ($self, $tree) = @_;
    my $leaves = $tree->get_all_leaves();
    for my $leaf (@$leaves) {
        mya_to_dtp_1path($leaf);
    }
}

sub mya_to_dtp_1path {
    my ($node) = @_;
    my $d = 0;
    for (;;) {
        my $dtp = 0;
        if ($node->has_tag('revised')) {
            if ($node->has_parent()) {
                $node = $node->parent();
                next;
            } else {
                return;
            }
        }
        if ($node->distance_to_parent != 0 && $node->has_parent) {
            $dtp = $node->distance_to_parent - $d;
        }
        $node->distance_to_parent($dtp);
        $node->add_tag("revised", "1");
        $d += $dtp;
        if ($node->has_parent()) {
            $node = $node->parent();
        } else {
            return;
        }
    }
}


sub ultrametrize {
    my ($self, $tree) = @_;
    my $longest_path = get_longest_path($tree);
    my $leaves = $tree->get_all_leaves();
    for my $leaf (@$leaves) {
        my $path = path_length($leaf);
        $leaf->distance_to_parent($leaf->distance_to_parent() + ($longest_path-$path));
    }
}

sub get_longest_path {
    my ($tree) = @_;
    my $leaves = $tree->get_all_leaves();
    my @paths;
    my $longest = -1;
    for my $leaf(@$leaves) {
        my $newpath = path_length($leaf);
        if ($newpath > $longest) {
            $longest = $newpath;
        }
    }
    return $longest;
}

sub binarize {
    my ($self, $orig_tree) = @_;
    my $newTree = $orig_tree->new();
    $newTree->name($orig_tree->name());
    $newTree->taxon_id($orig_tree->taxon_id);
    $newTree->genome_db_id($orig_tree->genome_db_id);
    $newTree->node_id($orig_tree->node_id());
    _binarize($orig_tree, $newTree, {});
    return $newTree;
}

sub _binarize {
    my ($origTree, $binTree, $taxon_ids) = @_;
    my $children = $origTree->children();
    for my $child (@$children) {
        my $newNode = $child->new();
        $newNode->taxon_id($child->taxon_id);
        $taxon_ids->{$child->parent->taxon_id}++;
        $newNode->genome_db_id($child->genome_db_id);
        $newNode->node_id($child->node_id());
        $newNode->distance_to_parent($child->distance_to_parent()); # no parent!!
        $newNode->node_name($child->name);
        if (scalar @{$binTree->children()} > 1) {
            $child->disavow_parent();
            my $newBranch = $child->new();
            for my $c (@{$binTree->children()}) {
                $c->distance_to_parent(0);
                $newBranch->add_child($c);
            }
            $binTree->add_child($newBranch);
            $newBranch->taxon_id($newBranch->parent->taxon_id);
            $newBranch->genome_db_id($newBranch->parent->genome_db_id);
            my $suffix = $taxon_ids->{$newBranch->parent->taxon_id} ? ".dup" . $taxon_ids->{$newBranch->parent->taxon_id} : "";
            $newBranch->node_name($newBranch->parent->name . $suffix);
        }
        $binTree->add_child($newNode);
        _binarize($child, $newNode, $taxon_ids);
    }
}

sub fix_zeros {
    my ($self, $tree) = @_;
    my $leaves = $tree->get_all_leaves();
    for my $leaf (@$leaves) {
        fix_zeros_1($leaf);
    }
}

sub fix_zeros_1 {
    my ($node) = @_;
    my $to_add = 0;
    for (;;) {
        return unless ($node->has_parent());
        my $dtp = $node->distance_to_parent();
        if ($dtp == 0) {
            $to_add++;
            $node->distance_to_parent(1);
        }
        my $siblings = $node->siblings;
        die "too many siblings" if (scalar @$siblings > 1);
        $siblings->[0]->distance_to_parent($siblings->[0]->distance_to_parent() + $to_add);
        $node = $node->parent();
    }
}

sub prune_tree {
    my ($self, $tree, $species_to_keep) = @_;

    my @nodes_to_remove = grep {!$species_to_keep->{$_->genome_db_id}} @{$tree->get_all_leaves};
    return $tree->remove_nodes(\@nodes_to_remove);
}


sub check_tree {
  my ($self, $tree) = @_;
  if (is_ultrametric($tree)) {
      if ($self->debug()) {
          print STDERR "The tree is ultrametric\n";
      }
  } else {
      die "The tree is NOT ultrametric\n";
  }

  is_binary($tree);
  if ($self->debug()) {
    print STDERR "The tree is binary\n";
  }
}

sub is_binary {
  my ($node) = @_;
  if ($node->is_leaf()) {
    return 0
  }
  my $children = $node->children();
  if (scalar @$children != 2) {
    my $name = $node->name();
    die "Not binary in node $name\n";
  }
  for my $child (@$children) {
    is_binary($child);
  }
}

sub is_ultrametric {
  my ($tree) = @_;
  my $leaves = $tree->get_all_leaves();
  my $path = -1;
  for my $leaf (@$leaves) {
    my $newpath = path_length($leaf);
    if ($path == -1) {
      $path = $newpath;
      next;
    }
    if ($path == $newpath) {
      $path = $newpath;
    } else {
      return 0
    }
  }
  return 1
}

sub path_length {
  my ($node) = @_;
  print STDERR "PATH LENGTH FOR ", $node->taxon_id;
  my $d = 0;
  for (;;){
    $d += $node->distance_to_parent();
    if ($node->has_parent()) {
      $node = $node->parent();
    } else {
      last;
    }
  }
  print STDERR " IS $d\n";
  return $d;
}

1;
