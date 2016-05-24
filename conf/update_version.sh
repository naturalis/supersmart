#!/bin/bash

# bash script to update version numbers as per 
# issue #41 before committing and kicking off a packer run

# shorter copy of SUPERSMART_HOME
H=$SUPERSMART_HOME

# semantic version tag as argument, MUST match /v\d+\.\d+\.\d+/ and MUST be at least
# one minor increment above the current highest git tag version.
TAGVERSION=$1
if [ -z "$1" ]; then 
	echo "Need varsion tag as argument. Must match /v\d+\.\d+\.\d+/"
	exit 0;
fi
eval $(perl -Mversion -e "version::is_strict('$TAGVERSION') or warn and print 'exit'");

# version in framework code
APIVERSION=$(perl -I$H/lib -MBio::SUPERSMART -e 'print $Bio::SUPERSMART::VERSION')

if [ "$TAGVERSION" != "$APIVERSION" ]; then

	# this so that we see what's happening
	echo "changing from version $APIVERSION to $TAGVERSION"

	# replace in the base module file
	echo "editing $H/lib/Bio/SUPERSMART.pm in place"
	sed -i '' -e "s/$APIVERSION/$TAGVERSION/" $H/lib/Bio/SUPERSMART.pm
	
	# replace in packer template
	echo "editing $H/conf/template.json in place"
	APIVERSION=$(echo $APIVERSION | sed -e 's/v//')
	TAGVERSION=$(echo $TAGVERSION | sed -e 's/v//')
	sed -i '' -e "s/$APIVERSION/$TAGVERSION/" $H/conf/template.json
	
	# commit tags and file changes
	git tag -m "Release v$TAGVERSION" v$TAGVERSION
	git push --tags
	git commit -m "Release v$TAGVERSION" $H/lib/Bio/SUPERSMART.pm $H/conf/template.json
	git push
	packer push template.json
else
	echo "version unchanged: (tag) $TAGVERSION == (api) $APIVERSION"
fi
