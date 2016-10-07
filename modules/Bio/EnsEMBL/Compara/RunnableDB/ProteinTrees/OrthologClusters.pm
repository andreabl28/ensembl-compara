=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Compara::RunnableDB::ProteinTrees::OrthologClusters

=head1 DESCRIPTION

This is the RunnableDB makes clusters consisting of given species list  orthologues(connected components, not necessarily complete graphs).

example:

standaloneJob.pl Bio::EnsEMBL::Compara::RunnableDB::ProteinTrees::OrthologClusters

=head1 AUTHORSHIP

Ensembl Team. Individual contributions can be found in the GIT log.

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with an underscore (_)

=cut

package Bio::EnsEMBL::Compara::RunnableDB::ProteinTrees::OrthologClusters;

use strict;
use warnings;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Compara::Graph::ConnectedComponentGraphs;
use base ('Bio::EnsEMBL::Compara::RunnableDB::GeneTrees::StoreClusters');
use Data::Dumper;

=head2
sub param_defaults {

    return {
            'species_set_id'         => 10000002, #[31,150,60], # macaque, human, orangutan
            'compara_db'           =>  'mysql://ensro@compara4/wa2_protein_trees_84',
            'reuse_db'				=> 'mysql://ensro@compara2/mp14_protein_trees_85',
            'sort_clusters'         => 1,  #needed by the store_clusterset sub 
            'member_type'           => 'protein', #needed by the store_clusterset sub 
            'mlss_id'				=> 40101, #needed by the store_clusterset sub


    };
}

=cut

=head2 fetch_input

	Description: pull orthologs for all pairwise combination of species in the list of given species

=cut

sub fetch_input {

	my $self = shift;
	$self->debug(5);
	my $species_set_id = $self->param('species_set_id');
#	print scalar $species_set_id;

	$self->param('previous_dba' , Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->go_figure_compara_dba($self->param('reuse_db')) );
	$self->param('prev_homolog_adaptor', $self->param('previous_dba')->get_HomologyAdaptor);
	$self->param('prev_ss_adaptor', $self->param('previous_dba')->get_SpeciesSetAdaptor);
	$self->param('prev_mlss_adaptor', $self->param('previous_dba')->get_MethodLinkSpeciesSetAdaptor);
	$self->param('mlss_adaptor', $self->compara_dba->get_MethodLinkSpeciesSetAdaptor);
	$self->param('homolog_adaptor', $self->compara_dba->get_HomologyAdaptor);
	$self->param('prev_ss_obj', $self->param('prev_ss_adaptor')->fetch_by_dbID($species_set_id) );
	my @gdb_objs = @ {$self->param('prev_ss_obj')->genome_dbs() };
	my @allOrthologs;

	for (my $gb1_index =0; $gb1_index < scalar @gdb_objs; $gb1_index++) {

		for (my $gb2_index = $gb1_index +1; $gb2_index < scalar @gdb_objs; $gb2_index++ ) {

			my $mlss = $self->param('prev_mlss_adaptor')->fetch_by_method_link_type_GenomeDBs('ENSEMBL_ORTHOLOGUES',
                       [ $gdb_objs[$gb1_index], $gdb_objs[$gb2_index] ]);
			print "\n   ", $mlss->dbID, "    mlss id \n" if $self->debug() ;
			my $homologs = $self->param('prev_homolog_adaptor')->fetch_all_by_MethodLinkSpeciesSet($mlss);
			print scalar @{ $homologs}, " homolog size \n" if $self->debug() ;
			push (@allOrthologs, @{$homologs});
			print scalar @allOrthologs, "  all size \n" if $self->debug() ;
		}
	} 
	$self->param('connected_split_genes', new Bio::EnsEMBL::Compara::Graph::ConnectedComponentGraphs);
	$self->param('ortholog_objects', \@allOrthologs );
}

sub run {
	my $self = shift;

    $self->_buildConnectedComponents($self->param('ortholog_objects'));

}

sub write_output {
    my $self = shift @_;

    $self->store_clusterset('default', $self->param('allclusters'));
}

sub _buildConnectedComponents {
	my $self = shift;
	my ($ortholog_objects) = @_;
    $self->dbc and $self->dbc->disconnect_if_idle();
    my $c = 0;
    my %allclusters = ();
    $self->param('allclusters', \%allclusters);
    while ( my $ortholog = shift( @{ $ortholog_objects } ) ) {
		my $gene_members = $ortholog->get_all_Members();
		my $seq_mid1 = $gene_members->[0]->dbID;
		my $seq_mid2 = $gene_members->[1]->dbID;
		print "seq mem ids   :   $seq_mid1     :    $seq_mid2   \n " if $self->debug() ;
		$self->param('connected_split_genes')->add_connection($seq_mid1, $seq_mid2);
		$c++;
#		last if $c >= 30;
		
	}
	printf("built %d clusters\n", $self->param('connected_split_genes')->get_graph_count) if $self->debug() ;
    printf("has %d distinct components\n", $self->param('connected_split_genes')->get_component_count) if $self->debug() ;
	my $cluster_id=0;
	my $holding_node = $self->param('connected_split_genes')->holding_node;
	foreach my $link (@{$holding_node->links}) {
    	my $one_node = $link->get_neighbor($holding_node);
    	my $nodes = $one_node->all_nodes_in_graph;
    	my @seq_member_ids = map {$_->node_id} @$nodes;
    	print Dumper(\@seq_member_ids) if $self->debug() ;
    	print "\n seq_member_ids :     \n" if $self->debug() ;
    	$allclusters{$cluster_id} = { 'members' => \@seq_member_ids };
    	$cluster_id++;	    

	}
}

1; 
