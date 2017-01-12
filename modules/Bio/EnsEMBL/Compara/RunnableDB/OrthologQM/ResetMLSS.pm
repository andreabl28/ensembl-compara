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

=pod

=head1 NAME
	
	Bio::EnsEMBL::Compara::RunnableDB::OrthologQM::ResetMLSS

=head1 SYNOPSIS

	Removes all scores in the ortholog_quality table assosiated with the
	list of input MLSS

=head1 DESCRIPTION

	Inputs:
	aln_mlss_ids	arrayref of method_link_species_set IDs to be removed

=cut

package Bio::EnsEMBL::Compara::RunnableDB::OrthologQM::ResetMLSS;

use strict;
use warnings;
use Data::Dumper;

use base ('Bio::EnsEMBL::Compara::RunnableDB::BaseRunnable');

sub run {
	my $self = shift;

	my $mlss_ids = $self->param('aln_mlss_ids');

	if ( !defined $mlss_ids || scalar( @{ $mlss_ids} ) < 1 ) {
		$self->warning("No MLSSs reset");
		return;
	}

	my $sql = 'DELETE FROM ortholog_quality WHERE alignment_mlss = ?';
	my $sth = $self->db->dbc->prepare($sql);
	foreach my $id ( @{ $mlss_ids } ){
		$sth->execute($id);
	}
}

1;