<?xml version="1.0"?>
<!DOCTYPE stylesheet [
	
]>

<!-- Customization layer -->
<xsl:stylesheet
	xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
	mlns:fo="http://www.w3.org/1999/XSL/Format"
	version="1.0">
	
	<!-- Use 'chunk.xsl' in line below to chunk files. -->
	<xsl:import href="/usr/share/sgml/docbook/xsl-stylesheets-1.61.2-2/xhtml/docbook.xsl"/>

	<!-- CSS -->
	<xsl:param name="html.stylesheet" select="'spanlib.css'"/>
	<xsl:param name="para.propagates.style" select="1"/>
	<xsl:param name="generate.id.attributes" select="1"/>
	 
	<!-- Shade -->
	<xsl:param name="shade.verbatim" select="1"/>
	<xsl:attribute-set name="shade.verbatim.style">
		<xsl:attribute name="bgcolor"></xsl:attribute>
		<xsl:attribute name="width">100%</xsl:attribute>
		<xsl:attribute name="align">center</xsl:attribute>
	</xsl:attribute-set>

	<xsl:template name="user.footer.content">
		<div class="generated">This document was generated <?dbtimestamp format="Y-m-d H:M:S"?> using xml/xsl and perl.</div>
	</xsl:template>

</xsl:stylesheet>
