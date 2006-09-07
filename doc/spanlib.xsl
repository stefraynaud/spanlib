<?xml version="1.0"?>
<!DOCTYPE stylesheet [
]>

<!-- Customization layer -->
<xsl:stylesheet
	xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
	xmlns:fo="http://www.w3.org/1999/XSL/Format"
	xmlns:date="http://exslt.org/dates-and-times"
	exclude-result-prefixes="date"
	version="1.0">

	<!-- Imports -->
	<xsl:import href="http://docbook.sourceforge.net/release/xsl/current/xhtml/chunk.xsl" />
	<xsl:import href="common.xsl"/>

	<!-- Chunking -->
	<!--xsl:param name="chunk.tocs.and.lots" select="1"/-->
	<xsl:param name="chunk.section.depth" select="2"/>
	<xsl:param name="chunker.output.indent" select="'yes'"/>

	<!-- CSS -->
	<xsl:param name="html.stylesheet" select="'spanlib.css'"/>

</xsl:stylesheet>
