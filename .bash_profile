# .bash_profile

# Get the aliases and functions
if [ -f ~/.bashrc ]; then
	. ~/.bashrc
fi


if [ -z "$PERL5LIB" ]
then
	# If PERL5LIB wasn't previously defined, set it...
	PERL5LIB=~/perl5reinstall/lib
else
	# ... otherwise, extend it.
	PERL5LIB=$PERL5LIB:~/perl5reinstall/lib
fi

MANPATH=$MANPATH:~/perl5reinstall/man 

# User specific environment and startup programs

GLOBUS_BIN=$HOME/globus/bin
PATH=$PATH:$HOME/bin:/usr/local/mvapich/bin:/usr/local/mpich2/bin:$GLOBUS_BIN
TMPDIR=$HOME/scratch/tmp

LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/nv/hp10/jstern7/globus/lib64:
LD_LIBRARY_PATH=$HOME/lib:$LD_LIBRARY_PATH

GLOBUS_LOCATION=$HOME/globus
MYPROXY_SERVER=myproxy.teragrid.org
MYPROXY_SERVER_PORT=7514


export PATH LD_LIBRARY_PATH PERL5LIB MANPATH TMPDIR GLOBUS_LOCATION MYPROXY_SERVER MYPROXY_SERVER_PORT GLOBUS_BIN

$GLOBUS_LOCATION/etc/globus-user-env.sh
