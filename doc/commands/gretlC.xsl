<?xml version="1.0"?>
<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform" >

  <!-- Stylesheet for C output -->

<xsl:param name="hlp">cli</xsl:param>
<xsl:param name="lang" select="'C'"/>

<xsl:output method="text" encoding="utf-8"/>

<xsl:variable name="intl"
	      select="document('hlp_l10n.xml')/internationalization"/>

<xsl:template name="gettext">
  <xsl:param name="key"/>
  <xsl:variable name="itext"
    select="$intl/localization[@language=$lang]/gentext[@key=$key]/@text"/>
  <xsl:choose>
    <xsl:when test="$itext">
      <xsl:value-of select="$itext"/>
    </xsl:when>
    <xsl:otherwise>
      <xsl:message terminate="yes">
        <xsl:text>** Error: no phrase with key = '</xsl:text>
        <xsl:value-of select="$key"/>
        <xsl:text>' found for lang '</xsl:text>
        <xsl:value-of select="$lang"/>
        <xsl:text>'.</xsl:text>
      </xsl:message>
    </xsl:otherwise>
  </xsl:choose>  
</xsl:template>

<xsl:template match="description" />
<xsl:template match="examples" />

<xsl:template match="funcref">
  <xsl:apply-templates/> 
</xsl:template>

<xsl:template match="funclist">
  <xsl:if test="@ref='functions'">
    <xsl:apply-templates/>
  </xsl:if>
</xsl:template>

<xsl:template match="function">
  <xsl:call-template name="gettext">
    <xsl:with-param name="key" select="@output"/>
  </xsl:call-template>
  <xsl:text> </xsl:text>
  <xsl:value-of select="@name"/>
  <xsl:choose>
    <xsl:when test="fnargs">
      <xsl:text> </xsl:text>
      <xsl:apply-templates/>
      <xsl:call-template name="nl"/>
    </xsl:when>
    <xsl:otherwise>
      <xsl:text> ()</xsl:text>
      <xsl:call-template name="dnl"/>
    </xsl:otherwise>
  </xsl:choose>
</xsl:template>

<xsl:template match="fnargs">
  <xsl:text>(</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>)</xsl:text>
  <xsl:call-template name="nl"/>
</xsl:template>

<xsl:template name="argname">
  <xsl:variable name="aname" select="current()"/>  
  <xsl:choose>
    <xsl:when test="contains($aname, '&amp;')">
      <xsl:value-of select ="substring($aname,2)"/>
    </xsl:when>
    <xsl:otherwise>
      <xsl:value-of select ="$aname"/>
    </xsl:otherwise>
  </xsl:choose>
</xsl:template>

<xsl:template match="fnarg">
  <xsl:call-template name="gettext">
    <xsl:with-param name="key" select="@type"/>
  </xsl:call-template>
  <xsl:if test="(@conditional)">
    <xsl:text> conditional</xsl:text>
  </xsl:if>
  <xsl:if test="@type != 'varargs'">
    <xsl:text> </xsl:text>
  </xsl:if>
  <xsl:call-template name="argname"/>
  <xsl:if test="(@optional)">
    <xsl:choose>
      <xsl:when test="@type='matrix' or @type='bundle' or
		      @type='string' or @type='matrixref' or
		      @type='bundleref' or @type='strings'">
	<xsl:text>[null]</xsl:text>
      </xsl:when>
      <xsl:otherwise>
	<xsl:text>[opt]</xsl:text>
      </xsl:otherwise>
    </xsl:choose>
  </xsl:if>  
  <xsl:if test="following-sibling::*">
    <xsl:text>, </xsl:text>
  </xsl:if>
</xsl:template>

<xsl:template match="func">
  <xsl:apply-templates/>
</xsl:template>

<xsl:template name="nl">
  <xsl:text>&#10;</xsl:text>  
</xsl:template>

<xsl:template name="dnl">
  <xsl:text>&#10;&#10;</xsl:text>  
</xsl:template>

</xsl:stylesheet>
