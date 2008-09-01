#!/usr/bin/perl -w

# Read in the transfac classification page and create a graph from that.



#<a name="1.1"></a>1.1 <i>Class</i>: <a href="/cgi-bin/pub/databases/transfac/getTF.cgi?AC=C0008">Leucine zipper factors (bZIP)</a>.
#<ul>
#<a name="1.1.1"></a>1.1.1 <i>Family</i>: AP-1(-like) components
#<ul>
#<a name="1.1.1.1"></a>1.1.1.1 <i>Subfamily</i>: Jun
#<ul>
#<a name="1.1.1.1.1"></a>1.1.1.1.1 <a href="/cgi-bin/pub/databases/transfac/getTF.cgi?AC=T00902">XBP-1</a> (human).
#<br>
#<a name="1.1.1.1.2"></a>1.1.1.1.2 <a href="/cgi-bin/pub/databases/transfac/getTF.cgi?AC=T00893">v-Jun</a> (ASV).




foreach $line (<STDIN>) {
	chomp $line;

	# Ignore blank lines
	if ($line eq "") {
		next;
	}


	if ($line =~ m/^<a name=\"(.+?)\"(.*)/) {
		$index = $1;
		$theRest = $2;

		if ($theRest =~ m/<i>Class<\/i>:/) {
		}
	}
}

