=head1 LICENSE

Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Compara::RunnableDB::ComparaHMM::BlastFactoryUnannotated

=head1 SYNOPSIS


=head1 DESCRIPTION

Fetch sorted list of member_ids and create jobs for BlastAndParsePAF. 
Supported keys:

   'step' => <number>
       How many sequences to write into the blast query file. Default 100

=cut

package Bio::EnsEMBL::Compara::RunnableDB::ComparaHMM::BlastFactoryUnannotated;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Compara::RunnableDB::BaseRunnable');


sub param_defaults {
    return {
            'step'  => 100,
    };
}

sub fetch_input {
    my ($self) = @_;

    my $sth = $self->compara_dba->get_HMMAnnotAdaptor->fetch_all_genes_missing_annot();
    my @unannotated_member_ids = sort {$a <=> $b} (map {$_->[0]} @{$sth->fetchall_arrayref});
    $sth->finish;
    $self->param('unannotated_member_ids', \@unannotated_member_ids);

}



sub write_output {
    my $self = shift @_;

    my $step = $self->param('step');
    my @member_id_list = @{$self->param('unannotated_member_ids')};

    while (@member_id_list) {
        my @job_array = splice(@member_id_list, 0, $step);
        my $output_id = {'start_member_id' => $job_array[0], 'end_member_id' => $job_array[-1]};
        $self->dataflow_output_id($output_id, 2);
    }
}

1;