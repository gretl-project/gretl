<?xml version="1.0"?>
<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform" >

<xsl:template match="commandlist"> 
 <xsl:apply-templates/> 
</xsl:template>

<xsl:template match="command">
<xsl:text>&#xa;#&#xa;</xsl:text>
<xsl:value-of select="@name"/>
<xsl:text>&#xa;@</xsl:text>
<xsl:value-of select="@section"/>
<xsl:text>&#xa;</xsl:text>
<xsl:apply-templates/>
<xsl:text>&#xa;&#xa;</xsl:text>
</xsl:template>

<xsl:template match="usage">
<xsl:apply-templates/>
<xsl:text>&#xa;</xsl:text>
</xsl:template>

<xsl:template match="options">
<xsl:text>&#xa;</xsl:text>
<xsl:apply-templates/>
</xsl:template>

<xsl:template match="option">
<xsl:choose>
  <xsl:when test="position() = 1"><xsl:text>&#xa;</xsl:text>
Options:    </xsl:when>
  <xsl:otherwise><xsl:text>            </xsl:text> 
  </xsl:otherwise>
</xsl:choose>
<xsl:apply-templates/>
<xsl:text>&#xa;</xsl:text>
</xsl:template>

<xsl:template match="arguments">
<xsl:text>&#xa;</xsl:text>
Arguments:  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="argument">
<xsl:if test="(@separated)">; </xsl:if> 
<xsl:apply-templates/><xsl:text> </xsl:text>
</xsl:template>

<xsl:template match="examples">
<xsl:text>&#xa;</xsl:text>
<xsl:apply-templates/>
</xsl:template>

<xsl:template match="example">
<xsl:choose>
  <xsl:when test="position() = 1">
Examples:   </xsl:when>
  <xsl:otherwise><xsl:text>            </xsl:text></xsl:otherwise>
</xsl:choose>
<xsl:apply-templates/>
<xsl:text>&#xa;</xsl:text>
</xsl:template>

<xsl:template match="flag">
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="effect">
  <xsl:text> (</xsl:text><xsl:apply-templates/><xsl:text>)</xsl:text>
</xsl:template>

<xsl:template match="repl">
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="cmd">
  <xsl:text>"</xsl:text><xsl:apply-templates/><xsl:text>"</xsl:text>
</xsl:template>

<xsl:template match="program">
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="literal">
  <xsl:text>"</xsl:text><xsl:apply-templates/><xsl:text>"</xsl:text>
</xsl:template>

<xsl:template match="mathvar">
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="book">
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="quote">
  <xsl:text>"</xsl:text><xsl:apply-templates/><xsl:text>"</xsl:text>
</xsl:template>

<xsl:template match="filename">
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="function">
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="equation">
<xsl:if test="(@status='display')">
<xsl:text>[/PARA]&#xa;&#xa;&#xa;  </xsl:text>
</xsl:if>
<xsl:value-of select="@ascii"/>
<xsl:if test="(@status='display')">
<xsl:text>&#xa;&#xa;&#xa;[PARA]</xsl:text>
</xsl:if>
</xsl:template>

<xsl:template match="para">
  <xsl:text>&#xa;[PARA]</xsl:text>
  <xsl:apply-templates/>[/PARA]
</xsl:template>

<xsl:template match="code">
  <xsl:text>&#xa;&#xa;</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>&#xa;&#xa;</xsl:text>
</xsl:template>

<xsl:template match="cmdref">
<xsl:text>"</xsl:text>
  <xsl:value-of select="@targ"/>
<xsl:text>"</xsl:text>
</xsl:template>

<xsl:template match="manref">
<xsl:value-of select="@before"/>
<xsl:text> the gretl manual </xsl:text>
<xsl:value-of select="@after"/>
</xsl:template>

<xsl:template match="menu-path">
Menu path:   <xsl:apply-templates/>
</xsl:template>

<xsl:template match="other-access">
Other access: <xsl:apply-templates/>
</xsl:template>

</xsl:stylesheet>
