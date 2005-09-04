<xsl:stylesheet version="1.0"
   xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

<xsl:import href="/usr/share/sgml/docbook-xsl/htmlhelp/htmlhelp.xsl"/>

<xsl:param name="htmlhelp.encoding" select="'UTF-8'"/>

<xsl:param name="lang" select="'en'"/>

<xsl:template match="command">
  <xsl:call-template name="inline.monoseq"/>
</xsl:template>

<xsl:template match="mathvar">
  <emphasis><xsl:apply-templates/></emphasis>
</xsl:template>

<xsl:param name="htmlhelp.hhp" select="'gretl.hhp'"/>
<xsl:param name="htmlhelp.chm" select="'gretl.chm'"/>
<xsl:param name="admon.graphics" select="1"/>
<xsl:param name="admon.graphics.extension" select="'.gif'"/>
<xsl:param name="admon.graphics.path">figures/</xsl:param>
<xsl:param name="navig.graphics" select="1"/>
<xsl:param name="graphic.default.extension" select="'gif'"/>

</xsl:stylesheet>


