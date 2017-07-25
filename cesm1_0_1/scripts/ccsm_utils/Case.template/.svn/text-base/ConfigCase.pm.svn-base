package ConfigCase;
my $pkg_nm = 'ConfigCase';
#-----------------------------------------------------------------------------------------------
#
# SYNOPSIS
# 
#   use ConfigCase;
# 
#   # read the configuration definition xml file
#   my $cfg = ConfigCase->new("$cfgdir/ccsm_utils/Case.template/config_definition.xml");
#
#   # set some parameters
#   $cfg->set($id, $value);
# 
#   # get some parameters
#   my $value = $cfg->get($id);
#
#   # Write an xml file out
#   $cfg->write_file("$caseroot/env_conf.xml", "xml");
#
#   # Write out documentation in a readme file
#   $cfg->write_doc("$caseroot/README/readme_env");
#
#   # Reset the config definition file with all of the values from the xml files
#   # in the $caseroot directory
#   $cfg->reset_setup("$caseroot/env_build.xml");
# 
# DESCRIPTION
# 
# ConfigCase objects are used to represent features of a CCSM model
# configuration that must be specified when a new case is created.
#
# new() Reads xml files that contain the configuration definition and
#       default values of the configuration parameters.
#
#       The "config_definition.xml" file contains all the allowable
#       parameters with a description of each one.  Where appropriate a
#       list of valid values of a parameter is given.  Generic default
#       values that are useful across many model configurations are
#       provided for some parameters.
#
#       The generic default values that are provided in the definition file
#       may be overridden by values in a setup file ("config_setup.xml")
#       that is assumed to provide appropriate values for a specific model
#       configuration.  This setup file is optional.
#
# is_valid_name() Returns true if the specified parameter name is contained in
#       the configuration definition file.
#
# is_ignore_name() Returns true if the specified parameter name is a name to ignore.
#
# is_valid_value() Returns true if the specified parameter name is contained in
#       the configuration definition file, and either 1) the specified value is
#       listed as a valid_value in the definition file, or 2) the definition file
#       doesn't specify the valid values.
#
# set() Sets values of the configuration parameters.  It takes the
#       parameter name and its value as arguments.  An invalid parameter
#       name (i.e., a name not present in the definition file) triggers an
#       exception.  If the definition file contains valid values of a
#       parameter, then the set method checks for a valid input value.  If
#       and invalid value is found then an exception is thrown.
#       ***NOTE*** If you don't want to trap exceptions, then use the query
#                  functions before calling this routine.
#
# get() Return the value of the specified configuration parameter.  Triggers
#       an exception if the parameter name is not valid.
#       ***NOTE*** If you don't want to trap exceptions, then use the query
#                  functions before calling this routine.
#
# write_file() Write a configuration xml file.  The first argument is the
#       filename.  The second argument, if present, is the commandline of the
#       configure command that was invoked to produce the output configuration
#       file.  It is written to the output file to help document the procedure
#       used to configure the executable.
#
# write_doc() Write documentation on the configuration to an output README file.
#
# reset_setup() Reset with all of the values from the xml files in the $caseroot directory
# 
# COLLABORATORS
# 
# IO::File
# XML::Lite
#-----------------------------------------------------------------------------------------------
#
# Date        Author                  Modification
#-----------------------------------------------------------------------------------------------
#  2008-Aug   Mariana Vertensten      Original version
#-----------------------------------------------------------------------------------------------

use strict;
use English;
#use warnings;
#use diagnostics;

use IO::File;
use XML::Lite;

sub new
{
    my $class = shift;
    my ($definition_file, $default_file) = @_;

    # bless the object here so the initialization has access to object methods
    my $cfg = {};
    bless( $cfg, $class );

    # Initialize the object with the configuration definition and its initial setup.
    $cfg->_initialize($definition_file, $default_file);

    return $cfg;
}

#-----------------------------------------------------------------------------------------------

sub is_ignore_name
{
# Return true if the requested name is a name to ignore
# These are descriptive names put in the config files -- but NOT in config_definition.xml

    my ($self, $name) = @_;

    if ( $name eq "NAME" or $name eq "SHORTNAME"    or $name eq "DESC" 
    or $name eq "STATUS" or $name eq "CCSM_PECOUNT" or $name eq "SUPPORT"
    or $name eq "GRID_MATCH" or $name eq "COMMENT"  or $name eq "SCIENCE" ) {
        return( 1 );
    } else {
        return( 0 );
    }
}

#-----------------------------------------------------------------------------------------------

sub is_valid_name
{
# Return true if the requested name is contained in the configuration definition.

    my ($self, $name) = @_;

    return defined($self->{$name}) ? 1 : 0;
}

#-----------------------------------------------------------------------------------------------

sub get
{
# Return requested value.

    my ($self, $name) = @_;

    defined($self->{$name}) or die "ERROR: unknown parameter name: $name\n";

    return $self->{$name}->{'value'};
}

#-----------------------------------------------------------------------------------------------

sub is_valid_value
{
# Return true if the specified parameter name is contained in
# the configuration definition file, and either 1) the specified value is
# listed as a valid_value in the definition file, or 2) the definition file
# doesn't specify the valid values.

    my ($self, $id, $value) = @_;

    # Check that the parameter name is in the configuration definition
    unless ($self->is_valid_name($id)) { return 0; }

    # Check that a list value is not supplied when parameter takes a scalar value.
    my $is_list_value = $self->{$id}->{'list'};
    unless ($is_list_value) {  # this conditional is satisfied when the list attribute is false, i.e., for scalars
	if ($value =~ /.*,.*/) { return 0; }   # the pattern matches when $value contains a comma, i.e., is a list
    }

    # Check that the value is valid
    my $valid_values = $self->{$id}->{'valid_values'};
    if ( $valid_values ne "" ) {  # if no valid values are specified, then $value is automatically valid
	if ($is_list_value) {
	    unless (_list_value_ok($value, $valid_values)) { return 0; }
	}
	else {
	    unless (_value_ok($value, $valid_values)) { return 0; }
	}

    }

    return 1;
}

#-----------------------------------------------------------------------------------------------

sub set
{
# Set requested value.
#
# This routine handles errors by throwing exceptions.  It will report exactly what problem was
# found in either the parameter name or requested value.
#
# To avoid dealing with exceptions use the is_valid_name(), is_valid_value() methods to get a
# true/false return before calling the set method.

    my ($self, $id, $value) = @_;

    # Check that the parameter name is in the configuration definition
    $self->is_valid_name($id) or die
	"ERROR: parameter name $id is not in the configuration definition\n";

    # Check that the value is valid
    my $valid_values = $self->{$id}->{'valid_values'};
    #$self->is_valid_value($id, $value) or die
	    #"ERROR: $value is not a valid value for parameter $id: valid values are $valid_values\n";

    # Get the type description hash for the variable and check that the type is valid
    # This method throws an exception when an error is encountered.
    my %type_ref = $self->_get_typedesc($id);
    $self->validate_variable_value($id, $value, \%type_ref);

    # Check that the value is valid
    if ( $valid_values ne "" ) {
       $self->is_valid_value($id, $value) or die
	       "ERROR: $value is not a valid value for parameter $id: valid values are $valid_values\n";
    }

    # Add the new value to the object's internal data structure.
    $self->{$id}->{'value'} = $value;

    return 1;
}

#-----------------------------------------------------------------------------------------------

sub write_file
{
# Write a configuration definition file.

    my $self = shift;
    my $filename = shift;   # filepath for output namelist
    my $format = shift;

    # determine what file to write
    my @groups;
    my @comps;
    my @models;
    my $xmode;

    $xmode = $self->get("XMLMODE");

    if ($filename =~ "env_case") {
	@groups = qw(case_def case_desc case_comp case_grid case_mach case_last);
    } elsif ($filename =~ "env_conf") {
	@groups = qw(conf_def conf_desc conf_map conf_ccsm conf_last); 
	@comps  = qw(cam datm clm dlnd cice dice pop2 docn cism);
	@models = qw(COMP_ATM COMP_LND COMP_ICE COMP_OCN COMP_GLC);
    } elsif ($filename =~ "env_mach_pes") {
	@groups = qw(mach_pes_def mach_pes_desc mach_pes_atm mach_pes_lnd mach_pes_ice mach_pes_ocn mach_pes_cpl mach_pes_glc mach_pes_stride mach_pes_last);
    } elsif ($filename =~ "env_build") {
	@groups = qw(build_def build_desc build_last);
    } elsif ($filename =~ "env_run") {
	@groups = qw(run_def run_desc run_stop run_rest run_flags run_diag run_cplhist run_cpl run_mach run_din run_dout run_last run_pop2);
       if ( $self->get("PTS_MODE") eq "TRUE" ) {
          push( @groups, "run_def_pts" );
       }
    }

    my $fh;
    if ( -f $filename ) { unlink( $filename ); }
    $fh = IO::File->new($filename, '>' ) or die "can't open file: $filename\n";

    if (($format eq "xml") || ($filename =~ m/xml/)) {
    print $fh <<"EOD";
<?xml version="1.0"?>

<config_definition>

EOD
    }

    if ($filename =~ "env_case.xml" && $xmode =~ "normal") {
    print $fh <<"EOD";
<!-- ========================================================================== -->
<!--                                                                            -->
<!--       These variables CANNOT BE CHANGED once a case has been created.      -->
<!--       Invoke create_newcase again if a different grid or component         -->
<!--       combination is required.                                             -->
<!--                                                                            -->
<!--       See README/readme_env and README/readme_general for details          -->
<!--                                                                            -->
<!-- ========================================================================== -->

EOD
}

    if ($filename =~ "env_conf.xml" && $xmode =~ "normal") {
    print $fh <<"EOD";
<!-- ========================================================================== -->
<!--                                                                            -->
<!--      These variables CANNOT be modified once configure -case has been      -->
<!--      invoked without first invoking configure -cleannamelist.              -->
<!--      These variables SHOULD NOT be modified once a run has been submitted. -->
<!--                                                                            -->
<!--      See README/readme_env and README/readme_general for details           -->
<!--                                                                            -->
<!-- ========================================================================== -->

EOD
}

    if ($filename =~ "env_run.xml" && $xmode =~ "normal") {
    print $fh <<"EOD";
<!-- ========================================================================== -->
<!--                                                                            -->
<!--       These variables MAY BE CHANGED ANYTIME during a run.                 -->
<!--       Additional machine speific variables that can be changed             -->
<!--       during a run are contained in the env_mach_specific file             -->
<!--                                                                            -->
<!--       Note1: users SHOULD NOT modify BUILD_COMPETE in env_build.xml        -->
<!--              this is done automatically by the scripts                     -->
<!--                                                                            -->
<!--       See README/readme_env and README/readme_general for details          -->
<!--                                                                            -->
<!-- ========================================================================== -->

EOD
}

    if ($filename =~ "env_build.xml" && $xmode =~ "normal") {
    print $fh <<"EOD";
<!-- ========================================================================== -->
<!--                                                                            -->
<!--       These variables SHOULD NOT be changed once the model has been built. -->
<!--       Currently, these variables are not cached.                           -->
<!--                                                                            -->
<!--       Note1: users SHOULD NOT modify BUILD_COMPETE below                   -->
<!--              this is done automatically by the scripts                     -->
<!--                                                                            -->
<!--       See README/readme_env and README/readme_general for details          -->
<!--                                                                            -->
<!-- ========================================================================== -->

EOD
}

    if ($filename =~ "env_mach_pes.xml" && $xmode =~ "normal") {
    print $fh <<"EOD";
<!-- ========================================================================== -->
<!--                                                                            -->
<!--      These variables CANNOT be modified once configure -case has been      -->
<!--      invoked without first invoking configure -cleanmach.                  -->
<!--                                                                            -->
<!--      See README/readme_env and README/readme_general for details           -->
<!--                                                                            -->
<!-- component task/thread settings                                             -->
<!-- if the user wants to change the values below after ./configure -case       -->
<!-- has been invoked, they can run                                             --> 
<!--    ./configure -cleanmach                                                  -->
<!--    ./configure -case                                                       -->
<!--  to reset the pes for the run                                              -->
<!--                                                                            -->
<!--  NTASKS are the total number of MPI tasks                                  -->
<!--  NTHRDS are the number of OpenMP threads per MPI task                      -->  
<!--  ROOTPE is the global mpi task associated with the root task               -->
<!--         of that component                                                  -->     
<!--  PSTRID is the stride of MPI tasks across the global                       -->
<!--         set of pes (for now this is set to 1)                              -->
<!--                                                                            -->
<!--  for example, for a setting with                                           -->  
<!--    NTASKS = 8                                                              -->
<!--    NTHRDS = 2                                                              -->
<!--    ROOTPE = 32                                                             -->
<!--  the MPI tasks would be placed starting on global pe 32                    -->
<!--  and each pe would be threaded 2-ways for this component.                  -->  
<!--                                                                            -->
<!--  Note: PEs that support threading never have an MPI task associated        -->
<!--    with them for performance reasons.  As a result, NTASKS and ROOTPE      -->
<!--    are relatively independent of NTHRDS and they determine                 -->
<!--    the layout of mpi processors between components.  NTHRDS is used        -->
<!--    to determine how those mpi tasks should be placed across the machine.   -->
<!--                                                                            -->
<!--  The following values should not be set by the user since they'll be       --> 
<!--  overwritten by scripts.                                                  -->
<!--    TOTALPES                                                                -->
<!--    CCSM_PCOST                                                              -->
<!--    CCSM_ESTCOST                                                            -->
<!--    PES_LEVEL                                                               -->
<!--    MAX_TASKS_PER_NODE                                                      -->
<!--    PES_PER_NODE                                                            -->
<!--    CCSM_TCOST                                                              -->
<!--    CCSM_ESTCOST                                                            -->
<!--                                                                            -->
<!--  The user can copy env_mach_pes.xml from another run, but they'll need to  -->
<!--  reconfigure and rebuild.                                                  -->
<!--    ./configure -cleanmach                                                  -->
<!--    ./configure -case                                                       -->
<!--    ./CASE.MACHINE.build                                                    -->
<!--                                                                            -->
<!-- ========================================================================== -->

EOD
}
    foreach my $group (@groups) {
	if (($format eq "xml") || ($filename =~ m/xml/)) {
	    if ($group =~ m/mach_pes_/ || $xmode =~ "expert") {
	        $self->_write_xml2($fh, $group);
	    } else {
		$self->_write_xml($fh, $group);
		print $fh "\n<!-- ====================================== -->";
	    }
	} else {
	    $self->_write_env($fh, $group);
	}
    }
    if ($self->get("COMP_OCN") eq "pop2") {
	if ($filename =~ "env_build") {
	    if (($format eq "xml") || ($filename =~ m/xml/)) {
               if ($xmode =~ "expert") {
   		   $self->_write_xml2($fh, "build_pop2");
	       } else {
		   $self->_write_xml($fh, "build_pop2");
		   print $fh "\n<!-- ====================================== -->";
	       }
	    } else {
		$self->_write_env($fh, "build_pop2");
	    }
	}
	if ($filename =~ "env_mach_pes") {
	    if (($format eq "xml") || ($filename =~ m/xml/)) {
		$self->_write_xml2($fh, "mach_pes_pop2");
	    } else {
		$self->_write_env($fh, "mach_pes_pop2");
	    }
	}
    }
    if ($self->get("COMP_ICE") eq "cice") {
	if ($filename =~ "env_build") {
	    if (($format eq "xml") || ($filename =~ m/xml/)) {
               if ($xmode =~ "expert") {
 		   $self->_write_xml2($fh, "build_cice");
	       } else {
 		   $self->_write_xml($fh, "build_cice");
		   print $fh "\n<!-- ====================================== -->";
	       }
	    } else {
		$self->_write_env($fh, "build_cice");
	    }
	}
	if ($filename =~ "env_mach_pes") {
	    if (($format eq "xml") || ($filename =~ m/xml/)) {
		$self->_write_xml2($fh, "mach_pes_cice");
	    } else {
		$self->_write_env($fh, "mach_pes_cice");
	    }
	}
    }
    if (@models && @comps) {
        my $comp  = undef;
        my $model = undef;
	foreach $comp (@comps)  {
	    foreach $model (@models) {
		if ($self->get($model) eq $comp) {
		    if (($format eq "xml") || ($filename =~ m/xml/)) {
                        if ($xmode =~ "expert") {
		 	   $self->_write_xml2($fh, "conf_$comp");
		       } else {
		 	   $self->_write_xml($fh, "conf_$comp");
			   print $fh "\n<!-- ====================================== -->";
		       }
		    } else {
			$self->_write_env($fh, "conf_$comp");
		    }
		}
	    }
	}
    }
    if (($format eq "xml") || ($filename =~ m/xml/)) {
    print $fh <<"EOD";

</config_definition>
EOD
    }
}

#-----------------------------------------------------------------------------------------------

sub write_doc
{
# Write the documentation on the configuration to an output README file.

    my $self = shift;
    my $filename = shift;   # filepath for output namelist

    my $fh;
    if ( -f $filename ) { unlink( $filename ); }
    $fh = IO::File->new($filename, '>' ) or die "can't open file: $filename\n";

    my @ids = keys %$self;
    foreach my $id (sort @ids) {

	print $fh "name: $id \n";
        my $ldesc = $self->{$id}->{'ldesc'};
        if ( ! defined($ldesc) ) { $ldesc = ""; }
	print $fh "description: " . $ldesc . "\n"; 
	if ($self->{$id}->{'valid_values'}) {
	    print $fh "valid values: $self->{$id}->{'valid_values'} \n"; 
	}
	print $fh "default value: $self->{$id}->{'value'} \n"; 
	if ( $self->{$id}->{'group'} =~ "case" ) {
	    print $fh "file: env_case.xml \n";
	    print $fh "locked_stage: create_newcase \n"; 
	}
	if ( $self->{$id}->{'group'} =~ "conf" ) {
	    print $fh "file: env_conf.xml \n";
	    print $fh "locked_stage: configure -case \n"; 
	}
	if ( $self->{$id}->{'group'} =~ "run" ) {
	    print $fh "file: env_run.xml \n";
	    print $fh "locked_stage: none \n"; 
	}
	if ( $self->{$id}->{'group'} =~ "build" ) {
	    print $fh "file: env_build.xml \n";
	    print $fh "locked_stage: none \n"; 
	}
	print $fh " \n";
    }
}

#-----------------------------------------------------------------------------------------------

sub write_docbook_master
{
# Write the documentation on the configuration to an output README file.

    my $self = shift;
    my $filename = shift;   # filepath for output namelist

    my $fh;
    if ( -f $filename ) { unlink( $filename ); }
    $fh = IO::File->new($filename, '>' ) or die "can't open file: $filename\n";

    my $gid;

    if ($filename =~ "case") { 
        $gid = "case";
	print $fh "<table><title>env_case.xml variables</title>\n";
#	print $fh "<tgroup cols=\"2\">\n";
#	print $fh "<tbody>\n";
    } elsif($filename =~ "conf") {
        $gid = "conf";
	print $fh "<table><title>env_conf.xml variables</title>\n";
#	print $fh "<tgroup cols=\"2\">\n";
#	print $fh "<thead>\n";
#	print $fh "<row> \n";
#	print $fh "<entry>Name</entry>\n";
#	print $fh "<entry>Description and [Valid Values] (Default)</entry>\n";
#	print $fh "</row> \n";
#	print $fh "</thead>\n";
#	print $fh "<tbody>\n";
    } elsif($filename =~ "build") {
        $gid = "build";
	print $fh "<table><title>env_build.xml variables</title>\n";
#	print $fh "<tgroup cols=\"2\">\n";
#	print $fh "<thead>\n";
#	print $fh "<row> \n";
#	print $fh "<entry>Name</entry>\n";
#	print $fh "<entry>Description and [Valid Values] (Default)</entry>\n";
#	print $fh "</row> \n";
#	print $fh "</thead>\n";
#	print $fh "<tbody>\n";
    } elsif($filename =~ "mach_pes") {
        $gid = "mach_pes";
	print $fh "<table><title>env_mach_pes.xml variables</title>\n";
#	print $fh "<tgroup cols=\"2\">\n";
#	print $fh "<thead>\n";
#	print $fh "<row> \n";
#	print $fh "<entry>Name</entry>\n";
#	print $fh "<entry>Description </entry>\n";
#	print $fh "</row> \n";
#	print $fh "</thead>\n";
#	print $fh "<tbody>\n";
    } elsif($filename =~ "run") {
        $gid = "run";
	print $fh "<table><title>env_run.xml variables</title>\n";
    }
	print $fh "<tgroup cols=\"4\">\n";
	print $fh "<thead>\n";
	print $fh "<row> \n";
	print $fh "<entry>Name</entry>\n";
	print $fh "<entry>Type</entry>\n";
	print $fh "<entry>Default</entry>\n";
	print $fh "<entry>Description [Valid Values]</entry>\n";
	print $fh "</row> \n";
	print $fh "</thead>\n";
	print $fh "<tbody>\n";
    
    my @ids = keys %$self;
    foreach my $id (sort @ids) {

        my $ldesc = "";
        my $valid = "";
        my $value = "";
        my $type  = "";

        $ldesc = $self->{$id}->{'ldesc'};
  	$valid = $self->{$id}->{'valid_values'};
	$value = $self->{$id}->{'value'};
	$type  = $self->{$id}->{'type'};

        if ( ! defined($ldesc) ) { $ldesc = ""; }
        if ( ! defined($valid) ) { $valid = ""; }
        if ( ! defined($value) ) { $value = ""; }
        if ( ! defined($type) ) { $type = ""; }

#	if ($filename =~ "case") { 
#	    if ( $self->{$id}->{'group'} =~ "case" ) {
#		print $fh "<row> \n";
#		print $fh "<entry>$id</entry>\n";
#		print $fh "<entry>$ldesc</entry>\n";
#		print $fh "</row>\n";
#	    }	    
#	} elsif ($filename =~ "conf") { 
#	    if ( $self->{$id}->{'group'} =~ "conf" ) {
#		print $fh "<row> \n";
#		print $fh "<entry>$id</entry>\n";
#		print $fh "<entry>Description: $ldesc. \n";
#		if ($self->{$id}->{'valid_values'}) {
#		    print $fh " Valid Values: [$self->{$id}->{'valid_values'}] ($self->{$id}->{'value'})</entry>\n";
#		} else {
#		    if ($self->{$id}->{'value'}) {
#			print $fh "Default: $self->{$id}->{'value'} </entry>\n";
#		    } else {
#			print $fh "</entry>\n";
#		    }
#		}
#		print $fh "</row>\n";
#	    }	    
#	} elsif ($filename =~ "mach_pes") { 
#	    if ( $self->{$id}->{'group'} =~ "mach_pes" ) {
#		print $fh "<row> \n";
#		print $fh "<entry>$id</entry>\n";
#		print $fh "<entry>$ldesc</entry>\n";
#		print $fh "</row>\n";
#	    }	    
#	} elsif ($filename =~ "run") { 
	    if ( $self->{$id}->{'group'} =~ $gid ) {
		print $fh "<row> \n";
		print $fh "<entry>$id</entry>\n";
		print $fh "<entry>$type</entry>\n";
		print $fh "<entry>$value</entry>\n";
		if ($self->{$id}->{'valid_values'}) {
		    print $fh "<entry>$ldesc [$valid] </entry>\n";
		} else {
		    print $fh "<entry>$ldesc</entry>\n";
                }
		print $fh "</row>\n";
	    }	    
#	} elsif ($filename =~ "build") { 
#	    if ( $self->{$id}->{'group'} =~ "build" ) {
#		print $fh "<row> \n";
#		print $fh "<entry>$id</entry>\n";
#		print $fh "<entry>Description: $ldesc. \n";
#		if ($self->{$id}->{'valid_values'}) {
#		    print $fh " Valid Values: [$self->{$id}->{'valid_values'}] ($self->{$id}->{'value'})</entry>\n";
#		} else {
#		    if ($self->{$id}->{'value'}) {
#			print $fh "Default: $self->{$id}->{'value'} </entry>\n";
#		    } else {
#			print $fh "</entry>\n";
#		    }
#		}
#		print $fh "</row>\n";
#	    }	    
#	}
    }
#    if ($filename =~ "case") { 
#	print $fh "</tbody>\n";
#	print $fh "</tgroup>\n";
#	print $fh "</table>\n";
#    } elsif($filename =~ "conf") {
#	print $fh "</tbody>\n";
#	print $fh "</tgroup>\n";
#	print $fh "</table>\n";
#    } elsif($filename =~ "build") {
#	print $fh "</tbody>\n";
#	print $fh "</tgroup>\n";
#	print $fh "</table>\n";
#    } elsif($filename =~ "mach_pes") {
#	print $fh "</tbody>\n";
#	print $fh "</tgroup>\n";
#	print $fh "</table>\n";
#    } elsif($filename =~ "run") {
	print $fh "</tbody>\n";
	print $fh "</tgroup>\n";
	print $fh "</table>\n";
#    }
    
}

#-----------------------------------------------------------------------------------------------

sub reset_setup
{

# Reset the config object from the setup file

    my ($self, $setup_file, $setup_id) = @_;

    # Process the setup file
    my $xml = XML::Lite->new( $setup_file );
    my $root = $xml->root_element();

    # Each parameter is contained in an "entry" element.  Get all these elements.
    my @elements = ();
    @elements = $xml->elements_by_name('entry');
    
    # Loop over the elements...
    foreach my $e (@elements) {
	
	# and extract the attributes
	my %attributes = $e->get_attributes();
	
	# just get the parameter name and value
	my $id    = $attributes{'id'};
	my $value = $attributes{'value'};
	
	# set new value
	if (defined($setup_id)) {
	    if ($id ne $setup_id) {$self->set($id, $value);}
	} else {
	    $self->set($id, $value);
	}
    } # end processing setup file
}

#-----------------------------------------------------------------------------------------------
# Private methods
#-----------------------------------------------------------------------------------------------

sub _initialize
{
# Read the configuration definition file.  Create an anonymous hash with the following
# structure:
# { id1 => {type=>'ttt', value=>"xxx", list =>"XXX", valid_values=>"yyy", 
#           group=>'vvv', sdesc=>'www', ldesc=>'WWW',definition=>"zzz"},
#   id2 => {type=>'ttt', value=>"xxx", list =>"XXX", valid_values=>"yyy", 
#           group=>'vvv', sdesc=>'www', ldesc=>'WWW',definition=>"zzz"},
#   ...												     										
#   idn => {type=>'ttt', value=>"xxx", list =>"XXX", valid_values=>"yyy", 
#           group=>'vvv', sdesc=>'www', ldesc=>'WWW',definition=>"zzz"},
# }

    my ($self, $definition_file, $setup_file) = @_;

    # Process the definition file
    my $xml = XML::Lite->new( $definition_file );
    my $root = $xml->root_element();

    # Check for valid root node
    my $name = $root->get_name();
    $name eq "config_definition" or die
	"ERROR: $definition_file is not a configuration definition file\n";

    # Each parameter is contained in an "entry" element.  Get all these elements.
    my @elements = $xml->elements_by_name('entry');

    # Loop over the elements...
    my $index = 0;
    foreach my $e (@elements) {

        # and extract the attributes and element content.
	my %attributes = $e->get_attributes();
	my $content    = $e->get_content();
	# if present strip initial and final newlines from content
	$content =~ s/^\n{1}//;
	$content =~ s/\n{1}$//;

	# Look for the specific attributes that are contained in the configuration definition.
	my $id    = $attributes{'id'};
	my $value = $attributes{'value'};
        my $sdesc = $attributes{'sdesc'};
        my $ldesc = $attributes{'ldesc'};
        my $list  = $attributes{'list'};
	my $group = $attributes{'group'};
	my $type  = $attributes{'type'};
	my $valid_values = defined $attributes{'valid_values'} ? $attributes{'valid_values'} : "";

	my $i = $index++; 
	# Now add the attributes and content to the object's internal data structure.
	$self->{$id} = {'type'         => $type, 
			'value'        => $value, 
			'list'         => $list, 
			'valid_values' => $valid_values, 
			'sdesc'        => $sdesc, 
			'ldesc'        => $ldesc, 
                        'definition'   => $content, 
			'group'        => $group,
                        'index'        => $i};
    }

    # Process the setup file
    if (defined $setup_file) {
	my $xml = XML::Lite->new( $setup_file );
	my $root = $xml->root_element();

	# Check for valid root node
	my $name = $root->get_name();
	$name eq "config_definition" or die
	    "ERROR: $definition_file is not a configuration definition file\n";

	# Each parameter is contained in an "entry" element.  Get all these elements.
	@elements = ();
	@elements = $xml->elements_by_name('entry');

	# Loop over the elements...
	foreach my $e (@elements) {

	    # and extract the attributes
	    my %attributes = $e->get_attributes();

	    # just get the parameter name and value
	    my $id    = $attributes{'id'};
	    my $value = $attributes{'value'};

	    # set new value
	    $self->set($id, $value);
	}
    } # end processing setup file
}

#-----------------------------------------------------------------------------------------------

sub _write_env
{
# Write a single line for the input env variable.

    my $self = shift;
    my $fh = shift;   # filepath for output namelist
    my $group = shift;

    print $fh <<"EOD";

EOD

    # first determine all of the indices you will write out
    my %id_indices;  
    my @ids = keys %$self;
    foreach my $id (sort @ids) {
	if ($group eq $self->{$id}->{'group'}) {
	    my $i = $self->{$id}->{'index'};
	    $id_indices{$i} = $id;
	}
    }

    # add the entry elements
    my @notsorted = keys %id_indices;
    my @sorted = sort { $a <=> $b } @notsorted;
    foreach my $i (@sorted) {
	my $id = $id_indices{$i};
	my $var = $self->{$id}->{'value'};
	if ($var =~ /\s/ ) { 
	    $var = "\"$var\"";
	} elsif ($var =~ /apos/ ) { 
	    $var = "\"$var\"";
	}
        $var =~ s/\&apos\;/\'/g;
        $var =~ s/\&lt\;/\</g;
	if ($group eq $self->{$id}->{'group'}) {
	    my $comment;
	    if ($self->{$id}->{'valid_values'}) {
		$comment = $self->{$id}->{'valid_values'};
	    } else {
		$comment = $self->{$id}->{'sdesc'};
	    }
	    print $fh <<"EOD";
setenv $id  $var     # $comment
EOD
        }
    }
}

#-----------------------------------------------------------------------------------------------

sub _write_xml
{
# Write a single XML element out to a file

    my $self = shift;
    my $fh = shift;   # filepath for output namelist
    my $group = shift;

    # separate the groups with spaces
    print $fh <<"EOD";

EOD

    # first determine all of the indices you will write out
    my %id_indices;  
    my @ids = keys %$self;
    foreach my $id (sort @ids) {
	if ($group eq $self->{$id}->{'group'}) {
	    my $i = $self->{$id}->{'index'};
            $id_indices{$i} = $id;
	}
    }

    # add the entry elements
    my @notsorted = keys %id_indices;
    my @sorted = sort { $a <=> $b } @notsorted;
    foreach my $i (@sorted) {
	my $id = $id_indices{$i};
	if ($self->{$id}->{'valid_values'} ne '') {
	    print $fh <<"EOD";

<!--"$self->{$id}->{'sdesc'}, valid values: $self->{$id}->{'valid_values'} ($self->{$id}->{'type'}) " -->
<entry id="$id"   value="$self->{$id}->{'value'}"  />    
EOD
} else {
		print $fh <<"EOD";

<!--"$self->{$id}->{'sdesc'} ($self->{$id}->{'type'}) " -->
<entry id="$id"   value="$self->{$id}->{'value'}"  />    
EOD
}

    }
}


#-----------------------------------------------------------------------------------------------

sub _write_xml2
{
# Write a single XML element out to a file

    my $self = shift;
    my $fh = shift;   # filepath for output namelist
    my $group = shift;

    # separate the groups with spaces
    print $fh <<"EOD";

EOD

    # first determine all of the indices you will write out
    my %id_indices;  
    my @ids = keys %$self;
    foreach my $id (sort @ids) {
	if ($group eq $self->{$id}->{'group'}) {
	    my $i = $self->{$id}->{'index'};
            $id_indices{$i} = $id;
	}
    }

    # add the entry elements
    my @notsorted = keys %id_indices;
    my @sorted = sort { $a <=> $b } @notsorted;
    foreach my $i (@sorted) {
	my $id = $id_indices{$i};
	print $fh <<"EOD";
<entry id="$id"   value="$self->{$id}->{'value'}"  />    
EOD
    }
}


#-----------------------------------------------------------------------------------------------

sub _list_value_ok
{
# Check that all input values ($values_in may be a comma separated list)
# are contained in the comma separated list of valid values ($valid_values).
# Return 1 (true) if all input values are valid, Otherwise return 0 (false).

    my ($values_in, $valid_values) = @_;

    my @values = split /,/, $values_in;

    my $num_vals = scalar(@values);
    my $values_ok = 0;

    foreach my $value (@values) {

	if (_value_ok($value, $valid_values)) { ++$values_ok; }

    }

    ($num_vals == $values_ok) ? return 1 : return 0;
}

#-----------------------------------------------------------------------------------------------

sub _value_ok
{

# Check that the input value is contained in the comma separated list of
# valid values ($valid_values).  Return 1 (true) if input value is valid,
# Otherwise return 0 (false).

    my ($value, $valid_values) = @_;

    # If the valid value list is null, all values are valid.
    unless ($valid_values) { return 1; }

    my @expect = split /,/, $valid_values;

    $value =~ s/^\s+//;
    $value =~ s/\s+$//;
    foreach my $expect (@expect) {
	if ($value =~ /^$expect$/) { return 1; }
    }

    return 0;
}

#-----------------------------------------------------------------------------------------------

sub get_var_type
{
# Return 'type' attribute for requested variable

    my ($self, $name) = @_;

    return $self->{$name}->{'type'};
}

#-----------------------------------------------------------------------------------------------

sub get_valid_values
{
# Return list of valid_values as an array for requested variable
# To return without quotes use the 'noquotes'=>1 option.
    my ($self, $name, %opts) = @_;

    my $valid_values = $self->{$name}->{'valid_values'};
    my $type = $self->{$name}->{'type'};
    my @values = split( /,/, $valid_values );

    # if string type and NOT noquote option and have a list -- add quotes around values
    if ( ! defined($opts{'noquotes'}) || ! $opts{'noquotes'} ) {
       if ( $#values > -1 && ($type eq "char") ) {
          for( my $i=0; $i <= $#values; $i++ ) {
             $values[$i] = "'$values[$i]'";
          }
       }
    }
    return( @values );
}

#-----------------------------------------------------------------------------------------------

sub _get_type
{
# Return 'type' attribute for requested variable

    my ($self, $name) = @_;

    return $self->{$name}->{'type'};
}

#-----------------------------------------------------------------------------------------------

sub _get_typedesc
#
# Return hash of description of data type read in from the file:
# Hash keys are:
#      type           type description (char, logical, integer, or real) (string)
#      validValues    Reference to array of valid values                 (string)
#
{
    my ($self, $name) = @_;
    my $nm = "_get_typedesc";

    my %datatype;
    my $type_def = $self->_get_type($name);
    if ($type_def =~ /^(char|logical|integer|real)/ ) {
      $datatype{'type'} = $1;
    } else {
	die "ERROR: in $nm (package $pkg_nm): datatype $type_def is NOT valid\n";
    }
    my @valid_values = $self->get_valid_values( $name );
    $datatype{'validValues'}  = \@valid_values;
    return( %datatype );
}

#-----------------------------------------------------------------------------------------------
# Perl regular expressions to match Fortran namelist tokens.

# Variable names.
# % for derived types, () for arrays
my $varmatch = "[A-Za-z_]{1,31}[A-Za-z0-9_]{0,30}[(0-9)%a-zA-Z_]*";

# Integer data.
my $valint = "[+-]?[0-9]+";
my $valint_repeat = "${valint}\\*$valint";

# Logical data.
my $vallogical1 = "[Tt][Rr][Uu][Ee]";
my $vallogical2 = "[Ff][Aa][Ll][Ss][Ee]";
my $vallogical = "$vallogical1|$vallogical2";
my $vallogical_repeat = "${valint}\\*$vallogical1|${valint}\\*$vallogical2";

# Real data.
# "_" are for f90 precision specification
my $valreal1 = "[+-]?[0-9]*\\.[0-9]+[EedDqQ]?[0-9+-]*";
my $valreal2 = "[+-]?[0-9]+\\.?[EedDqQ]?[0-9+-]*";
my $valreal = "$valreal1|$valreal2";
my $valreal_repeat = "${valint}\\*$valreal1|${valint}\\*$valreal2";

# Match for all valid data-types: integer, real or logical
# note: valreal MUST come before valint in this string to prevent integer portion of real 
#       being separated from decimal portion
my $valall = "$vallogical|$valreal|$valint";

# Match for all valid data-types with repeater: integer, real, logical, or string data
# note: valrepeat MUST come before valall in this string to prevent integer repeat specifier 
#       being accepted alone as a value
my $valrepeat = "$vallogical_repeat|$valreal_repeat|$valint_repeat";

# Match for all valid data-types with or without numberic repeater at the lead
my $valmatch = "$valrepeat|$valall";

# Same as above when a match isn't required
my $nrvalmatch = $valmatch. "||";

#-----------------------------------------------------------------------------------------------

sub validate_variable_value
#
# Validate that a given value matches the expected input type definition
# Expected description of keys for the input type hash is:
#      type           type description (char, logical, integer, or real)       (string)
#      validValues    Reference to array of valid values                       (string)
#
{
  my ($self, $var, $value, $type_ref) = @_;
  my $nm = "validate_variable_value";

  # Ensure type hash has required variables
  if ( ref($type_ref) !~ /HASH/ ) {
      die "ERROR: in $nm : Input type is not a HASH reference.\n";
  }
  foreach my $item ( "type", "validValues" ) {
     if ( ! exists($$type_ref{$item}) ) {
         die "ERROR: in $nm: Variable name $item not defined in input type hash.\n";
     }
  }

  # If not string -- check that array size is smaller than definition
  my $str_len = 0;
  if ( $$type_ref{'type'} eq "char" ) {
      $str_len = 1;
  }
  my @values;
  if ( $str_len == 0 ) {
      @values = split( /,/, $value );
      # Now check that values are correct for the given type
      foreach my $i ( @values ) {
	  my $compare;
	  if (      $$type_ref{'type'} eq "logical" ) {
	      $compare = $vallogical;
	  } elsif ( $$type_ref{'type'} eq "integer" ) {
	      $compare = $valint;
	  } elsif ( $$type_ref{'type'} eq "real" ) {
	      $compare = $valreal;
	  } else {
	      die "ERROR: in $nm (package $pkg_nm): Type of variable name $var is " . 
		  "not a valid FORTRAN type (logical, integer, real, or char).\n";
	  }
	  if ( $i !~ /^\s*(${compare})$/ ) {
	      die "ERROR: in $nm (package $pkg_nm): Variable name $var " .
		  "has a value ($i) that is not a valid type " . $$type_ref{'type'} . "\n";
	  }
      }
  }
}

#-----------------------------------------------------------------------------------------------

1; # to make use or require happy
