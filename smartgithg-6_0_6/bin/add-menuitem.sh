#!/bin/bash
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

# absolutize dir
oldpwd=`pwd`
cd "${SMARTGIT_BIN}"
SMARTGIT_BIN=`pwd`
cd "${oldpwd}"

ICON_NAME=syntevo-smartgithg-6
TMP_DIR=`mktemp --directory`
DESKTOP_FILE=$TMP_DIR/syntevo-smartgithg-6.desktop
cat << EOF > $DESKTOP_FILE
[Desktop Entry]
Version=1.0
Encoding=UTF-8
Name=SmartGit/Hg 6
Keywords=git;hg;svn;mercurial;subversion
GenericName=Git&Hg-Client + SVN-support
Type=Application
Categories=Development;RevisionControl
Terminal=false
StartupNotify=true
Exec="$SMARTGIT_BIN/smartgithg.sh" %u
MimeType=x-scheme-handler/smartgit
Icon=$ICON_NAME.png
EOF

xdg-desktop-menu install $DESKTOP_FILE
xdg-icon-resource install --size  32 "$SMARTGIT_BIN/smartgithg-32.png"  $ICON_NAME
xdg-icon-resource install --size  48 "$SMARTGIT_BIN/smartgithg-48.png"  $ICON_NAME
xdg-icon-resource install --size  64 "$SMARTGIT_BIN/smartgithg-64.png"  $ICON_NAME
xdg-icon-resource install --size 128 "$SMARTGIT_BIN/smartgithg-128.png" $ICON_NAME

rm $DESKTOP_FILE
rm -R $TMP_DIR
