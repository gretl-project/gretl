<?xml version="1.0"?>
<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform" >

  <!-- Stylesheet for C output -->

<xsl:param name="hlp">cli</xsl:param>
<xsl:param name="lang" select="'en'"/>

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
  <xsl:text>"</xsl:text>
  <xsl:value-of select="@name"/>
  <xsl:text>", return: </xsl:text>
  <xsl:value-of select="@output"/>
  <xsl:text>, </xsl:text>
  <xsl:apply-templates/>
  <xsl:call-template name="nl"/>
</xsl:template>

<xsl:template match="fnargs">
  <xsl:text>argc=</xsl:text>
  <xsl:value-of select="count(fnarg)"/>
  <xsl:text>,</xsl:text>
  <xsl:call-template name="nl"/>
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="fnarg">
  <xsl:choose>
    <xsl:when test="@type='varargs' or @type='seebelow'">
      <xsl:text>varargs</xsl:text>
    </xsl:when>
    <xsl:otherwise>
      <xsl:apply-templates/>
    </xsl:otherwise>
  </xsl:choose>
  <xsl:text> (</xsl:text>
  <xsl:call-template name="gettext">
    <xsl:with-param name="key" select="@type"/>
  </xsl:call-template>
  <xsl:if test="(@optional)">
    <xsl:choose>
      <xsl:when test="@type='matrixref'">
        <xsl:text>,</xsl:text>
        <xsl:call-template name="gettext">
          <xsl:with-param name="key" select="'or'"/>
        </xsl:call-template>
        <xsl:text>null</xsl:text>
      </xsl:when>
      <xsl:otherwise>
        <xsl:text>, </xsl:text>
        <xsl:call-template name="gettext">
          <xsl:with-param name="key" select="'optional'"/>
        </xsl:call-template>
      </xsl:otherwise>
    </xsl:choose>
  </xsl:if>
  <xsl:if test="(@conditional)">
    <xsl:text>, </xsl:text>
    <xsl:call-template name="gettext">
      <xsl:with-param name="key" select="'conditional'"/>
    </xsl:call-template>
  </xsl:if> 
  <xsl:text>)</xsl:text>
  <xsl:text>&#10;</xsl:text>
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
