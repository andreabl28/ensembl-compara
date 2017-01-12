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

=head1 NAME

Bio::EnsEMBL::Compara::PipeConfig::EPO_pt2_conf

=head1 SYNOPSIS

    #1. Update ensembl-hive, ensembl and ensembl-compara GIT repositories before each new release

    #3. Check all default_options, you will probably need to change the following :
        pipeline_db (-host)
        resource_classes 

	'ensembl_cvs_root_dir' - the path to the compara/hive/ensembl GIT checkouts - set as an environment variable in your shell
        'password' - your mysql password
	'compara_anchor_db' - database containing the anchor sequences (entered in the anchor_sequence table)
	'compara_master' - location of your master db containing relevant info in the genome_db, dnafrag, species_set, method_link* tables
        The dummy values - you should not need to change these unless they clash with pre-existing values associated with the pairwise alignments you are going to use

    #4. Run init_pipeline.pl script:
        Using command line arguments:
        init_pipeline.pl Bio::EnsEMBL::Compara::PipeConfig::EPO_pt2_conf.pm

    #5. Run the "beekeeper.pl ... -sync" and then " -loop" command suggested by init_pipeline.pl

    #6. Fix the code when it crashes
=head1 DESCRIPTION  

    This configuaration file gives defaults for mapping (using exonerate at the moment) anchors to a set of target genomes (dumped text files)

=head1 CONTACT

Please email comments or questions to the public Ensembl
developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

Questions may also be sent to the Ensembl help desk at
<http://www.ensembl.org/Help/Contact>.

=cut

package Bio::EnsEMBL::Compara::PipeConfig::EPO_pt2_conf;

use strict;
use warnings;

use Bio::EnsEMBL::Hive::Version 2.4;

use base ('Bio::EnsEMBL::Compara::PipeConfig::ComparaGeneric_conf');

sub default_options {
    my ($self) = @_;

    return {
	%{$self->SUPER::default_options},

        # Change this name
        'species_set_name' => '17mammals_reuse',

        # But not this one
        'pipeline_name' => $self->o('species_set_name').'_epo_anchor_mapping',

        # Where the pipeline lives
        'host'  => 'compara2',

	  # database containing the anchors for mapping
	'compara_anchor_db' => 'mysql://ensro@compara3/sf5_TEST_gen_anchors_mammals_cat_100',

	'mapping_exe' => "/software/ensembl/compara/exonerate/exonerate",
	'mapping_params'    => { bestn=>11, gappedextension=>"no", softmasktarget=>"no", percent=>75, showalignment=>"no", model=>"affine:local", },
	 # place to dump the genome sequences
	'seq_dump_loc' => '/data/blastdb/Ensembl/' . 'compara_genomes_test_' . $self->o('ensembl_release'),
	 # dont dump the MT sequence for mapping
	'dont_dump_MT' => 1,
	 # batch size of grouped anchors to map
	'anchor_batch_size' => 10,
	 # max number of sequences to allow in an anchor
	'anc_seq_count_cut_off' => 15,
        'compara_master' => 'mysql://ensro@compara1/mm14_ensembl_compara_master',
     };
}

sub pipeline_create_commands {
    my ($self) = @_; 
    return [
        @{$self->SUPER::pipeline_create_commands},  # inheriting database and hive tables' creation
        'mkdir -p '.$self->o('seq_dump_loc'),
           ];  
}

sub resource_classes {
    my ($self) = @_; 
    return {
	%{$self->SUPER::resource_classes},  # inherit 'default' from the parent class
	'default' => {'LSF' => '-C0 -M2500 -R"select[mem>2500] rusage[mem=2500]"' }, # farm3 lsf syntax
	'mem3500' => {'LSF' => '-C0 -M3500 -R"select[mem>3500] rusage[mem=3500]"' },
	'mem7500' => {'LSF' => '-C0 -M7500 -R"select[mem>7500] rusage[mem=7500]"' },
	'mem14000' => {'LSF' => '-C0 -M14000 -R"select[mem>14000] rusage[mem=14000]"' },

    };  
}

sub pipeline_wide_parameters {
	my $self = shift @_;
	return {
		%{$self->SUPER::pipeline_wide_parameters},

                'epo_mlss_id' => $self->o('epo_mlss_id'),
		'seq_dump_loc' => $self->o('seq_dump_loc'),
		'compara_anchor_db' => $self->o('compara_anchor_db'),
		'master_db' => $self->o('compara_master'),
		'reuse_db' => $self->o('reuse_db'),
	};
	
}

sub pipeline_analyses {
	my ($self) = @_;
	print "pipeline_analyses\n";

        return [

# ---------------------------------------------[copy tables from master and fix the offsets]---------------------------------------------

            {   -logic_name => 'copy_table_factory',
                -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
                -parameters => {
                    'db_conn'      => '#master_db#',
                    'inputlist'    => [ 'method_link', 'ncbi_taxa_name', 'ncbi_taxa_node' ],
                    'column_names' => [ 'table' ],
                },
                -input_ids => [{}],
                -flow_into => {
                    '2->A' => { 'copy_table' => { 'src_db_conn' => '#db_conn#', 'table' => '#table#' } },
                    'A->1' => [ 'load_genomedb_factory' ],
                },
            },

            {   -logic_name    => 'copy_table',
                -module        => 'Bio::EnsEMBL::Hive::RunnableDB::MySQLTransfer',
                -parameters    => {
                    'mode'          => 'topup',
                    'filter_cmd'    => 'sed "s/ENGINE=MyISAM/ENGINE=InnoDB/"',
                },
            },


# ---------------------------------------------[load GenomeDB entries from master+cores]---------------------------------------------

            {   -logic_name => 'load_genomedb_factory',
                -module     => 'Bio::EnsEMBL::Compara::RunnableDB::GenomeDBFactory',
                -parameters => {
                    'mlss_id'               => $self->o('mlss_id'),
                    'compara_db'            => '#master_db#',   # that's where genome_db_ids come from
                    'extra_parameters'      => [ 'locator' ],
                },
                -flow_into => {
                    '2->A' => { 'load_genomedb' => { 'master_dbID' => '#genome_db_id#', 'locator' => '#locator#' }, },
                    'A->1' => [ 'create_mlss_ss' ],
                },
            },

            {   -logic_name => 'load_genomedb',
                -module     => 'Bio::EnsEMBL::Compara::RunnableDB::LoadOneGenomeDB',
                -parameters => {
                    'registry_conf_file' => $self->o('reg_conf'),
                },
                -flow_into => [ 'copy_dnafrags_from_master' ],
            },

            {   -logic_name => 'copy_dnafrags_from_master',
                -module     => 'Bio::EnsEMBL::Hive::RunnableDB::MySQLTransfer',
                -parameters => {
                    'src_db_conn'   => '#master_db#',
                    'table'         => 'dnafrag',
                    'where'         => 'genome_db_id = #genome_db_id#',
                    'mode'          => 'insertignore',
                },
                -flow_into => [ 'check_reusability' ],
            },

            {   -logic_name => 'check_reusability',
                -module     => 'Bio::EnsEMBL::Compara::RunnableDB::ProteinTrees::CheckGenomedbReusability',
                -parameters => {
                    check_gene_content  => 0,
                },
                -flow_into => {
                    2 => '?accu_name=reused_gdb_ids&accu_address=[]&accu_input_variable=genome_db_id',
                    3 => '?accu_name=nonreused_gdb_ids&accu_address=[]&accu_input_variable=genome_db_id',
                },
            },

            {   -logic_name => 'create_mlss_ss',
                -module     => 'Bio::EnsEMBL::Compara::RunnableDB::GeneTrees::PrepareSpeciesSetsMLSS',
                -parameters => {
                    'tree_method_link' => 'EPO',
                    'create_homology_mlss'  => 0,
                },
                -flow_into => [ 'reuse_anchor_align_factory' ],
            },

            {   -logic_name     => 'reuse_anchor_align_factory',
                -module     => 'Bio::EnsEMBL::Compara::RunnableDB::GenomeDBFactory',
                -parameters => {
                    'species_set_id'    => '#reuse_ss_id#',
                },
                -flow_into => {
                    '2->A' => [ 'reuse_anchor_align' ],
                    'A->1' => [ 'dump_genome_sequence_factory' ],
                },
            },

            # Copy all the untrimmed anchor_aligns
            {   -logic_name => 'reuse_anchor_align',
                -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
                -parameters => {
                    'db_conn'    => '#reuse_db#',
                    # FIXME the INSERT may fail because method_link_species_set_id is not defined in this database
                    'inputquery' => 'SELECT anchor_align.* FROM anchor_align JOIN dnafrag USING (dnafrag_id) WHERE genome_db_id = #genome_db_id# AND untrimmed_anchor_align_id IS NULL',
                },
                -flow_into => {
                    2 => [ '?table_name=anchor_align' ],
                },
            },

            {   -logic_name     => 'dump_genome_sequence_factory',
                -module     => 'Bio::EnsEMBL::Compara::RunnableDB::GenomeDBFactory',
                -parameters => {
                    'species_set_id'    => '#nonreuse_ss_id#',
                    'extra_parameters'      => [ 'locator', 'name', 'assembly' ],
                },
                -flow_into => {
                    '2->A' => { 'dump_genome_sequence' => { 'genome_db_name' => '#name#', 'genome_db_assembly' => '#assembly#', 'genome_db_id' => '#genome_db_id#' } },
                    'A->1'  => [ 'remove_overlaps' ],
                },
            },

	    {	-logic_name     => 'dump_genome_sequence',
		-module         => 'Bio::EnsEMBL::Compara::Production::EPOanchors::DumpGenomeSequence',
		-parameters => {
			'dont_dump_MT' => $self->o('dont_dump_MT'),
		},
		-flow_into => [ 'map_anchors_factory' ],
		-rc_name => 'mem7500',
		-hive_capacity => 10,
	    },

	    {	-logic_name     => 'map_anchors_factory',
                -module         => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
		-parameters     => {
                    'db_conn'       => '#compara_anchor_db#',
                    'inputquery'    => 'SELECT anchor_id FROM anchor_sequence GROUP BY anchor_id HAVING COUNT(*) <= #anc_seq_count_cut_off# ORDER BY anchor_id',
                    'contiguous'    => 0,
                    'step'          => $self->o('anchor_batch_size'),
                    'anc_seq_count_cut_off' => $self->o('anc_seq_count_cut_off'),
		},
		-flow_into => {
			2 => { 'map_anchors' => { 'anchor_ids' => '#_range_list#', 'genome_db_file' => '#genome_db_file#', 'genome_db_id' => '#genome_db_id#' } },
		},
		-rc_name => 'mem7500',
		-hive_capacity => 10,
	    },

	    {	-logic_name     => 'map_anchors',
		-module         => 'Bio::EnsEMBL::Compara::Production::EPOanchors::MapAnchors',
		-parameters => {
			'mapping_exe' => $self->o('mapping_exe'),
	                'mapping_params' => $self->o('mapping_params'),
		},
		-hive_capacity => 1000,
		-max_retry_count => 1,
	    },

	    {	-logic_name     => 'remove_overlaps',
		-module         => 'Bio::EnsEMBL::Compara::Production::EPOanchors::RemoveAnchorOverlaps',
		-rc_name => 'mem3500',
		-flow_into => {
			1 => [ 'trim_anchor_align_factory' ],
		},
	    },

            {   -logic_name => 'trim_anchor_align_factory',
                -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
                -parameters => {
                                'inputquery'      => "SELECT DISTINCT(anchor_id) AS anchor_id FROM anchor_align WHERE untrimmed_anchor_align_id IS NULL AND is_overlapping = 0",
                               },  
                -flow_into => {
                               2 => [ 'trim_anchor_align' ],
                              },  
		-rc_name => 'mem3500',
            },  

	    {   -logic_name => 'trim_anchor_align',			
		-module     => 'Bio::EnsEMBL::Compara::Production::EPOanchors::TrimAnchorAlign',
		-parameters => {
				'method_link_species_set_id' => '#epo_mlss_id#',
			},
                -flow_into => {
                    -1 => 'trim_anchor_align_himem',
                },
		-hive_capacity => 150,
		-batch_size    => 20,
	    },

	    {   -logic_name => 'trim_anchor_align_himem',
		-module     => 'Bio::EnsEMBL::Compara::Production::EPOanchors::TrimAnchorAlign',
		-parameters => {
				'method_link_species_set_id' => '#epo_mlss_id#',
			},
		-hive_capacity => 150,
                -rc_name => 'mem3500',
	    },
    ];
}	


1;
