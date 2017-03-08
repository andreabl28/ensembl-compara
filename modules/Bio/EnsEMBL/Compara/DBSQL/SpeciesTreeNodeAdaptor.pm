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

Bio::EnsEMBL::Compara::DBSQL::SpeciesTreeNodeAdaptor

=head1 DESCRIPTION

  SpeciesTreeNodeAdaptor - Adaptor for different species trees used in ensembl-compara


=head1 APPENDIX

  The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Compara::DBSQL::SpeciesTreeNodeAdaptor;

use strict;
use warnings;
use Data::Dumper;


use Bio::EnsEMBL::Compara::SpeciesTreeNode;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Scalar qw(:assert);

use base ('Bio::EnsEMBL::Compara::DBSQL::NestedSetAdaptor', 'Bio::EnsEMBL::Compara::DBSQL::TagAdaptor');

=head2 new_from_NestedSet

    Arg[1]      : An object that inherits from NestedSet
    Example     : my $st_node = Bio::EnsEMBL::Compara::SpeciesTreeNode->new_from_NestedSet($tree);
    Description : Constructor for species tree nodes. Given an object that inherits from NestedSet (possibly a tree), creates a new SpeciesTreeNode (possibly a tree).
    ReturnType  : EnsEMBL::Compara::SpeciesTreeNode
    Exceptions  : none
    Caller      : General

=cut

sub new_from_NestedSet {
    my ($self, $nestedSet_tree) = @_;

    my $genomeDB_Adaptor = $self->db->get_GenomeDBAdaptor;
    my $NCBITaxon_Adaptor = $self->db->get_NCBITaxonAdaptor;

    my $tree = $nestedSet_tree->cast('Bio::EnsEMBL::Compara::SpeciesTreeNode');
    for my $node (@{$tree->get_all_nodes}) {
        if ($node->is_leaf && !$node->genome_db_id) {
            my $genomeDB;
            if ($node->genome_db_id) {
                $genomeDB = $genomeDB_Adaptor->fetch_by_dbID($node->genome_db_id);
            } elsif (defined $node->taxon_id) {
                $genomeDB = $genomeDB_Adaptor->fetch_by_taxon_id($node->taxon_id);
            } elsif (defined $node->node_name) {
                $genomeDB = $genomeDB_Adaptor->fetch_by_name_assembly($node->node_name);
            }
            if (defined $genomeDB) {
                $node->genome_db_id($genomeDB->dbID);
            }
        }
        if (!$node->{_taxon} || !$node->{_taxon_id}) {
            my $taxon_node;
            if (defined $node->taxon_id) {
                $taxon_node = $NCBITaxon_Adaptor->fetch_node_by_taxon_id($node->taxon_id);
            } elsif (defined $node->name) {
                $taxon_node = $NCBITaxon_Adaptor->fetch_node_by_name($node->name);
            }
            if (defined $taxon_node) {
                $node->taxon($taxon_node);
            }
        }
    }
    return $tree;
}


## TODO: This is very similar to GeneTreeNodeAdaptor's store_node, maybe we can
## abstract out this code in NestedSetAdaptor
sub store_node {
    my ($self, $node, $mlss_id) = @_;

## This may fail in the case of CAFEGeneFamilyNodes
    assert_ref($node, 'Bio::EnsEMBL::Compara::SpeciesTreeNode', 'node');

    if (not($node->adaptor and $node->adaptor->isa('Bio::EnsEMBL::Compara::DBSQL::SpeciesTreeNodeAdaptor') and $node->adaptor eq $self)) {

        my $count_current_rows = $self->dbc->db_handle->selectall_arrayref("SELECT COUNT(*) FROM species_tree_node")->[0]->[0];

        my $node_id;
        if ($count_current_rows) {
            $node_id = $self->generic_insert('species_tree_node', {}, 'node_id');

        } else {
            ## if table is empty then $node_id should be mlss_id * 1000
            $node_id = $mlss_id * 1000;
            $self->generic_insert('species_tree_node', {'node_id' => $node_id});
        }

        $node->node_id( $node_id );
    }

    my $rc = $self->generic_update('species_tree_node',
        {
            'parent_id'             => $node->parent ? $node->parent->node_id : undef,
            'root_id'               => $node->root->node_id,
            'left_index'            => $node->left_index,
            'right_index'           => $node->right_index,
            'distance_to_parent'    => $node->distance_to_parent,
            'taxon_id'              => $node->taxon_id,
            'genome_db_id'          => $node->genome_db_id,
            'node_name'             => $node->node_name,
        }, {
            'node_id'               => $node->node_id,
        } );
    if ($rc == 0) {
        die "Could not update the newly-created species_tree_node row node_id=".$node->node_id.". Please investigate\n";
    }

    return $node->node_id;
}


#
# tagging
#
sub _tag_capabilities {
    return ('species_tree_node_tag', 'species_tree_node_attr', 'node_id', 'node_id', 'tag', 'value');
}


#################################################
#
# subclass override methods
#
#################################################

sub _columns {
    return qw ( stn.node_id
                stn.root_id
                stn.parent_id
                stn.left_index
                stn.right_index
                stn.distance_to_parent
                stn.taxon_id
                stn.genome_db_id
                stn.node_name
             );
}

sub _tables {
    return (['species_tree_node', 'stn'], ['species_tree_root','str']);
}

sub _default_where_clause {
    return "stn.root_id = str.root_id";
}

sub create_instance_from_rowhash {
    my ($self, $rowhash) = @_;

    my $node = new Bio::EnsEMBL::Compara::SpeciesTreeNode;

    $self->init_instance_from_rowhash($node, $rowhash);
    return $node;
}

sub init_instance_from_rowhash {
    my ($self, $node, $rowhash) = @_;

    # SUPER is NestedSetAdaptor
    $self->SUPER::init_instance_from_rowhash($node, $rowhash);

    $node->taxon_id($rowhash->{'taxon_id'});
    $node->genome_db_id($rowhash->{'genome_db_id'});
    $node->node_name($rowhash->{'node_name'});

    $node->adaptor($self);
    return $node;
}

1;
