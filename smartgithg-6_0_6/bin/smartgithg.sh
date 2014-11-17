#!/bin/bash
#
# Normally, editing this script should not be required.
#
# To specify an alternative Java Runtime Environment, set the environment variable SMARTGITHG_JAVA_HOME

case "$BASH" in
    */bash) :
        ;;
    *)
        exec /bin/bash "$0" "$@"
        ;;
esac

if [ "$SMARTGITHG_JAVA_HOME" = "" ] ; then
	SMARTGITHG_JAVA_HOME=$SMARTGIT_JAVA_HOME
fi
if [ "$SMARTGITHG_JAVA_HOME" = "" ] && [ -f "/usr/lib/jvm/java-7-openjdk-i386/jre/bin/java" ] ; then
	SMARTGITHG_JAVA_HOME="/usr/lib/jvm/java-7-openjdk-i386/jre"
fi
if [ "$SMARTGITHG_JAVA_HOME" = "" ] ; then
	SMARTGITHG_JAVA_HOME=$JAVA_HOME
fi

if [ "$SMARTGITHG_MAX_HEAP_SIZE" = "" ] ; then
	SMARTGITHG_MAX_HEAP_SIZE=$SMARTGIT_MAX_HEAP_SIZE
fi
if [ "$SMARTGITHG_MAX_HEAP_SIZE" = "" ] ; then
	SMARTGITHG_MAX_HEAP_SIZE=256m
fi

# this seems necessary for Solaris to find the Cairo-library
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/usr/lib/gnome-private/lib

_JAVA_EXEC="java"
if [ "$SMARTGITHG_JAVA_HOME" != "" ] ; then
    _TMP="$SMARTGITHG_JAVA_HOME/bin/java"
    if [ -f "$_TMP" ] ; then
        if [ -x "$_TMP" ] ; then
            _JAVA_EXEC="$_TMP"
        else
            echo "Warning: $_TMP is not executable"
        fi
    else
        echo "Warning: $_TMP does not exist"
    fi
fi

if ! which "$_JAVA_EXEC" >/dev/null ; then
    echo "Error: No java environment found"
    exit 1
fi

#
# Resolve the location of the SmartGit/Hg installation.
# This includes resolving any symlinks.
PRG=$0
while [ -h "$PRG" ]; do
    ls=`ls -ld "$PRG"`
    link=`expr "$ls" : '^.*-> \(.*\)$' 2>/dev/null`
    if expr "$link" : '^/' 2> /dev/null >/dev/null; then
        PRG="$link"
    else
        PRG="`dirname "$PRG"`/$link"
    fi
done

SMARTGIT_BIN=`dirname "$PRG"`

# Without the following line sliders are not visible in Ubuntu 12.04
# (see <https://bugs.eclipse.org/bugs/show_bug.cgi?id=368929>)
export LIBOVERLAY_SCROLLBAR=0

# absolutize dir
oldpwd=`pwd`
cd "${SMARTGIT_BIN}"; SMARTGIT_BIN=`pwd`
cd "${oldpwd}"; unset oldpwd

SMARTGIT_HOME=`dirname "$SMARTGIT_BIN"`


# work-around for SWT bug https://bugs.eclipse.org/bugs/show_bug.cgi?id=419729
if type "lsb_release" > /dev/null 2> /dev/null ; then
    UBUNTU_VERSION=`lsb_release -a 2>/dev/null | grep Release | cut -d ':' -f2 | tr -d '\t' | tr -d ' '`

    if [ "$UBUNTU_VERSION" == "13.10" ] || [ "$UBUNTU_VERSION" == "14.04" ] || [ "$UBUNTU_VERSION" == "14.10" ] ; then
        UBUNTU_MENUPROXY=0
    fi
fi


_VM_PROPERTIES="-Dsun.io.useCanonCaches=false"

# Uncomment the following line to change the location where SmartGit/Hg should store
# settings (the given example path will make SmartGit/Hg portable by storing the settings
# in the installation directory):
#_VM_PROPERTIES="$_VM_PROPERTIES -Dsmartgit.settings=\${smartgit.installation}/.smartgit"

$_JAVA_EXEC $_VM_PROPERTIES -Xmx${SMARTGITHG_MAX_HEAP_SIZE} -Xverify:none -Dsmartgit.vm-xmx=${SMARTGITHG_MAX_HEAP_SIZE} -jar "$SMARTGIT_HOME/lib/bootloader.jar" "$@"
