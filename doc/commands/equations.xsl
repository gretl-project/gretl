<?xml version="1.0"?>
<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform" >

<xsl:output method="text" encoding="iso-8859-1"/>

<xsl:template match="equation">
  <xsl:text>&lt;texequation </xsl:text>
  <xsl:text>fileref="</xsl:text>
  <xsl:value-of select="@graphic"/>
  <xsl:text>"&gt;&#10;</xsl:text>
  <xsl:value-of select="normalize-space(@tex)"/>
  <xsl:text>&#10;&lt;/texequation&gt;&#10;</xsl:text>
</xsl:template>

<xsl:template match="commandlist">
  <xsl:text>&lt;equation-set latexopt="12pt" density="96x96" usepackage="mathptmx"&gt;&#10;</xsl:text>
  <xsl:apply-templates select=".//equation"/>
  <xsl:text>&lt;/equation-set&gt;&#10;</xsl:text>
</xsl:template>

</xsl:stylesheet>
