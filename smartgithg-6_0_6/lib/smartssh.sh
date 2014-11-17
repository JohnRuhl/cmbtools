#!/bin/sh
"${SMARTGIT_JAVA_HOME}/bin/java" -cp "${SMARTGIT_CLASSPATH}" -Dsmartgit.logging=false -Djava.net.preferIPv4Stack=true com.syntevo.dvcs.transport.ssh.SdSshMain "$@"
exit 0
