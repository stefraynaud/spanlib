<?xml version="1.0"?>
<!DOCTYPE stylesheet [
<!ENTITY css SYSTEM "spanlib.css">
]>

<!-- Customization layer -->
<xsl:stylesheet
	xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
	xmlns:fo="http://www.w3.org/1999/XSL/Format"
	xmlns:date="http://exslt.org/dates-and-times"
	exclude-result-prefixes="date"
	version="1.0">

	<!-- Imports -->
	<xsl:import href="http://docbook.sourceforge.net/release/xsl/current/xhtml/docbook.xsl" />
	<xsl:import href="common.xsl"/>

	<!-- Out formating -->
	<xsl:output indent="yes" />

	<!-- CSS -->
	<xsl:template name="user.head.content">
		<style type="text/css">
			&css;
		</style>
	</xsl:template>


</xsl:stylesheet>
