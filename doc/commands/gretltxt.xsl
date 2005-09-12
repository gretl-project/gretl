<?xml version="1.0"?>
<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform" >

<xsl:param name="hlp">cli</xsl:param>
<xsl:param name="lang" select="'en'"/>

<xsl:output method="text" encoding="iso-8859-1"/>

<xsl:variable name="phrases"
  select="document('hlpstrings.xml')/phrases"/>

<xsl:template name="gettext">
  <xsl:param name="key"/>
  <xsl:variable name="phrase"
    select="$phrases/phrase[@key=$key and @lang=$lang]"/>
  <xsl:choose>
    <xsl:when test="$phrase">
      <xsl:value-of select="$phrase"/>
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

<xsl:template match="commandlist">
  <xsl:if test="$hlp='gui'">
    <xsl:text>@new-style gretl.hlp&#10;</xsl:text>     
  </xsl:if>
<xsl:apply-templates/> 
</xsl:template>

<xsl:template match="command">
  <xsl:if test="not(@context) or @context=$hlp">
  <xsl:if test="position() > 1">
    <xsl:call-template name="nl"/>
  </xsl:if>
  <xsl:choose>
    <xsl:when test="$hlp='gui'">
      <xsl:text># </xsl:text>
      <xsl:value-of select="@name"/>
      <xsl:text> </xsl:text>
      <xsl:value-of select="@section"/>
      <xsl:text> "</xsl:text>
      <xsl:value-of select="@label"/>
      <xsl:text>"&#10;</xsl:text>
    </xsl:when>
    <xsl:otherwise>
      <xsl:text>#&#10;</xsl:text>
      <xsl:value-of select="@name"/>
      <xsl:text>&#10;@</xsl:text>
      <xsl:value-of select="@section"/>
      <xsl:if test="not(usage)">
        <xsl:call-template name="nl"/>
      </xsl:if>
    </xsl:otherwise>
  </xsl:choose>
  <xsl:apply-templates/>
  <xsl:call-template name="dnl"/>
  <xsl:if test="(not(@context) and $hlp='gui')">
    <xsl:call-template name="gettext">
      <xsl:with-param name="key" select="'scriptcommand'"/>
    </xsl:call-template>
    <xsl:value-of select="@name"></xsl:value-of>
    <xsl:call-template name="nl"/>
  </xsl:if>
</xsl:if>
</xsl:template>

<xsl:template match="description">
  <xsl:if test="not(@context) or @context=$hlp">
    <xsl:apply-templates/>
  </xsl:if>
</xsl:template>

<xsl:template match="usage">
  <xsl:if test="$hlp='cli'">
    <xsl:apply-templates/>
    <xsl:call-template name="dnl"/>
  </xsl:if>
</xsl:template>

<xsl:template match="arguments">
  <xsl:text>&#xa;</xsl:text>
  <xsl:choose>
    <xsl:when test="count(argument) > 1">
      <xsl:call-template name="gettext">
        <xsl:with-param name="key" select="'args'"/>
      </xsl:call-template>
    </xsl:when>
    <xsl:otherwise>
      <xsl:call-template name="gettext">
        <xsl:with-param name="key" select="'arg'"/>
      </xsl:call-template>
    </xsl:otherwise> 
  </xsl:choose>
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="argblock">
  <xsl:if test="(@separated)">; </xsl:if>
  <xsl:if test="(@optional)">[ </xsl:if>
  <xsl:apply-templates/>
  <xsl:if test="(@separated)">; </xsl:if>
  <xsl:if test="(@optional)">] </xsl:if>
</xsl:template>

<xsl:template match="argument">
  <xsl:if test="(@optional)">[ </xsl:if> 
  <xsl:if test="(@separated) and not(preceding-sibling::*/@separated)">; </xsl:if>
  <xsl:if test="(@alternate)">
    <xsl:call-template name="gettext">
      <xsl:with-param name="key" select="'or'"/>
    </xsl:call-template>
  </xsl:if>
  <xsl:if test="@flag">
    <xsl:value-of select="@flag"/>
  </xsl:if> 
  <xsl:apply-templates/>
  <xsl:text> </xsl:text>
  <xsl:if test="(@optional)">] </xsl:if> 
</xsl:template>

<xsl:template match="options">
  <xsl:call-template name="nl"/>
  <xsl:choose>
    <xsl:when test="count(option) > 1">
      <xsl:call-template name="gettext">
        <xsl:with-param name="key" select="'opts'"/>
      </xsl:call-template>
    </xsl:when>
    <xsl:otherwise>
      <xsl:call-template name="gettext">
        <xsl:with-param name="key" select="'opt'"/>
      </xsl:call-template>
    </xsl:otherwise> 
  </xsl:choose>
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="option|example">
  <xsl:if test="position() > 1">
    <xsl:text>&#xa;            </xsl:text>
  </xsl:if> 
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="altforms">
  <xsl:call-template name="nl"/>
  <xsl:call-template name="gettext">
    <xsl:with-param name="key" select="'altforms'"/>
  </xsl:call-template>
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="altform">
  <xsl:if test="position() > 1">
    <xsl:text>&#xa;                  </xsl:text>
  </xsl:if> 
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="option|example">
  <xsl:if test="position() > 1">
    <xsl:text>&#xa;            </xsl:text>
  </xsl:if> 
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="examples">
  <xsl:call-template name="nl"/>
  <xsl:choose>
    <xsl:when test="count(example) > 1">
      <xsl:call-template name="gettext">
        <xsl:with-param name="key" select="'examples'"/>
      </xsl:call-template>
    </xsl:when>
    <xsl:otherwise>
      <xsl:call-template name="gettext">
        <xsl:with-param name="key" select="'example'"/>
      </xsl:call-template>
    </xsl:otherwise> 
  </xsl:choose>
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="flag">
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="effect">
  <xsl:text> (</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>)</xsl:text>
</xsl:template>

<xsl:template match="repl">
  <xsl:if test="@quote='true'">
    <xsl:text>"</xsl:text>
  </xsl:if>
  <xsl:apply-templates/>
  <xsl:if test="@quote='true'">
    <xsl:text>"</xsl:text>
  </xsl:if>
</xsl:template>

<xsl:template match="cmd">
  <xsl:text>"</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>"</xsl:text>
</xsl:template>

<xsl:template match="program">
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="lit">
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="mathvar">
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="book">
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="quote">
  <xsl:text>"</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>"</xsl:text>
</xsl:template>

<xsl:template match="filename">
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="func">
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
  <xsl:if test="not(@context) or @context=$hlp">
  <xsl:choose>
    <xsl:when test="parent::li and ancestor::ilist">
      <xsl:text>&#xa;[ILISTPAR]</xsl:text>
      <xsl:apply-templates/>[/ILISTPAR]
    </xsl:when>
    <xsl:when test="parent::li and ancestor::nlist">
      <xsl:text>&#xa;[NLISTPAR]</xsl:text>
      <xsl:if test="count(../*)=1">
        <xsl:value-of select="1 + count(../preceding-sibling::*)"/>
        <xsl:text>. </xsl:text>
      </xsl:if>
      <xsl:apply-templates/>[/NLISTPAR]
    </xsl:when>
    <xsl:otherwise>
      <xsl:text>&#xa;[PARA]</xsl:text>
      <xsl:apply-templates/>[/PARA]
    </xsl:otherwise>
  </xsl:choose> 
 </xsl:if>
</xsl:template>

<xsl:template match="code">
  <xsl:if test="not(@context) or @context=$hlp">
    <xsl:call-template name="dnl"/>
    <xsl:apply-templates/>
    <xsl:call-template name="dnl"/>
  </xsl:if>
</xsl:template>

<xsl:template match="ilist">
  <xsl:if test="not(@context) or @context=$hlp">
    <xsl:apply-templates/>
  </xsl:if>
</xsl:template>

<xsl:template match="nlist">
  <xsl:if test="not(@context) or @context=$hlp">
    <xsl:apply-templates/>
  </xsl:if>
</xsl:template>

<xsl:template match="li">
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="cmdref">
<xsl:text>"</xsl:text>
  <xsl:value-of select="@targ"/>
<xsl:text>"</xsl:text>
</xsl:template>

<xsl:template match="manref">
  <xsl:call-template name="gettext">
    <xsl:with-param name="key" select="'fullman'"/>
  </xsl:call-template>
</xsl:template>

<xsl:template match="tabref">
  <xsl:call-template name="gettext">
    <xsl:with-param name="key" select="'tablebelow'"/>
  </xsl:call-template>
</xsl:template>

<xsl:template match="menu-path">
  <xsl:text>&#xa;</xsl:text>
  <xsl:call-template name="gettext">
    <xsl:with-param name="key" select="'menupath'"/>
  </xsl:call-template>
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="other-access">
  <xsl:text>&#xa;</xsl:text>
  <xsl:call-template name="gettext">
    <xsl:with-param name="key" select="'otheraccess'"/>
  </xsl:call-template>
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="table"> 
<xsl:if test="not(@context) or @context=$hlp">
  <xsl:text>[TABLE]&#xa;</xsl:text>
  <xsl:if test="@id">
    <xsl:text>[ROW]&#xa;[CELL]&#xa;</xsl:text>
    <xsl:value-of select="@lhead"/>
    <xsl:text>[/CELL]&#xa;[CELL]&#xa;</xsl:text>
    <xsl:value-of select="@rhead"/>
    <xsl:text>[/CELL]&#xa;[/ROW]&#xa;</xsl:text>
    <xsl:text>[ROW]&#xa;[CELL]&#xa;-------</xsl:text>
    <xsl:text>[/CELL]&#xa;[CELL]&#xa;-------</xsl:text>
    <xsl:text>[/CELL]&#xa;[/ROW]&#xa;</xsl:text>
  </xsl:if>
  <xsl:apply-templates/>
  <xsl:text>[/TABLE]&#xa;</xsl:text>  
</xsl:if>
</xsl:template>

<xsl:template match="row">
  <xsl:text>[ROW]&#xa;</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>[/ROW]&#xa;</xsl:text>
</xsl:template>

<xsl:template match="cell">
  <xsl:text>[CELL]&#xa;</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>[/CELL]&#xa;</xsl:text>
</xsl:template>

<xsl:template match="footnote"/>

<xsl:template name="nl">
  <xsl:text>&#10;</xsl:text>  
</xsl:template>

<xsl:template name="dnl">
  <xsl:text>&#10;&#10;</xsl:text>  
</xsl:template>

</xsl:stylesheet>
