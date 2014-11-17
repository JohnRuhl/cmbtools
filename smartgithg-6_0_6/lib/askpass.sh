#!/bin/sh
"${SMARTGIT_JAVA_HOME}/bin/java" -cp "${SMARTGIT_CLASSPATH}" -Dsmartgit.logging=false -Djava.net.preferIPv4Stack=true com.syntevo.smartgit.transport.askpass.SgAskPasswordMain "$@"
exit 0
