<?xml version="1.0"?>
<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform" >

  <xsl:output method="text" encoding="utf-8"/>
  <xsl:strip-space elements="*" />

<xsl:template match="commandref">
  <xsl:apply-templates/> 
</xsl:template>

<xsl:template match="usage" />
<xsl:template match="gui-access" />
<xsl:template match="subhead" />

<xsl:template match="command">
  <xsl:if test="@context='gui'">
    <xsl:apply-templates select="description"/>
  </xsl:if>
</xsl:template>

<xsl:template match="description">
  <xsl:apply-templates select="para"/>
</xsl:template>

<xsl:template match="para">
  <xsl:apply-templates select="cite"/>
</xsl:template>

<xsl:template match="cite">
  <xsl:text>&lt;@bib="</xsl:text>
  <xsl:value-of select="normalize-space()"/>
  <xsl:text>;</xsl:text>
  <xsl:value-of select="@key"/>
  <xsl:text>"&gt;</xsl:text>
  <xsl:text>&#10;</xsl:text>
</xsl:template>

</xsl:stylesheet>
