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

	<!-- Automatic numbering -->
	<xsl:param name="section.autolabel" select="1"/>
	<xsl:param name="section.autolabel.max.depth" select="4"/>
	<xsl:param name="appendix.autolabel" select="1"/>
	<xsl:param name="appendix.autolabel.max.depth" select="3"/>

	<!-- Misc -->
	<xsl:param name="use.id.as.filename" select="1"/>
	<xsl:param name="para.propagates.style" select="1"/>

	<!-- Shade -->
	<!--xsl:param name="shade.verbatim" select="1"/-->
	<xsl:attribute-set name="shade.verbatim.style">
		<xsl:attribute name="width">100%</xsl:attribute>
		<xsl:attribute name="align">center</xsl:attribute>
	</xsl:attribute-set>

	<!-- Additional footer -->
	<xsl:template name="user.footer.navigation">
		<a href="http://sourceforge.net"><img src="http://sflogo.sourceforge.net/sflogo.php?group_id=168272&amp;type=1" width="88" height="31" border="0" alt="SourceForge.net Logo" /></a>
		<div class="generated">This document was generated using xsltproc and perl.</div>
	</xsl:template>

	<!-- Date -->
<xsl:template match="pubdate[@role='now']" mode="titlepage.mode">
     <xsl:variable name="now" select="date:date-time()"/>
     <fo:block>
       <xsl:text>This documentation was generated on </xsl:text>
       <xsl:value-of select="date:day-in-month($now)"/>
       <xsl:text> </xsl:text>
       <xsl:value-of select="date:month-name($now)"/>
       <xsl:text> </xsl:text>
       <xsl:value-of select="date:year($now)"/>
       <!--xsl:text> at </xsl:text>
       <xsl:value-of select="date:hour-in-day($now)"/>
       <xsl:text>:</xsl:text>
       <xsl:value-of select="date:minute-in-hour($now)"/-->
     </fo:block>
   </xsl:template>

	<!-- Graphics be used for admonitions (notes, warnings) -->
	<!--xsl:param name="admon.graphics" select="1"/>
	<xsl:param name="admon.graphics.path">images/</xsl:param>
	<xsl:param name="admon.graphics.extension" select="'.gif'"/-->


	<!-- Date in header -->
	<!--xsl:template name="user.head.content">
	<meta name="date">
		<xsl:attribute name="content">
			<xsl:call-template name="datetime.format">
				<xsl:with-param name="date" select="date:date-time()"/>
				<xsl:with-param name="format" select="'c'"/>
			</xsl:call-template>
		</xsl:attribute>
	</meta>
	</xsl:template-->

</xsl:stylesheet>
