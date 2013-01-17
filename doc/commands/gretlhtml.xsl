<?xml version="1.0"?>
<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform" >

  <!-- Stylesheet for HTML output -->

<xsl:param name="hlp">cli</xsl:param>
<xsl:param name="topic">cmds</xsl:param>
<xsl:param name="lang" select="'en'"/>

<xsl:output method="text" encoding="utf-8"/>

<xsl:variable name="intl"
  select="document('hlp_l10n.xml')/internationalization"/>

<xsl:template name="gettext">
  <xsl:param name="key"/>
  <xsl:variable name="itext"
    select="concat(normalize-space($intl/localization[@language=$lang]/gentext[@key=$key]/@text),' ')"/>
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

<xsl:template name="gettext-nospace">
  <xsl:param name="key"/>
  <xsl:variable name="itext"
    select="normalize-space($intl/localization[@language=$lang]/gentext[@key=$key]/@text)"/>
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

<xsl:template name="index-cell">
  <xsl:param name="entry"/>
  <xsl:text>&lt;td&gt;</xsl:text>
  <xsl:text>&lt;a href="#</xsl:text>
  <xsl:value-of select="$entry"/>
  <xsl:text>"&gt;</xsl:text>
  <xsl:value-of select="$entry"/>
  <xsl:text>&lt;/a&gt;</xsl:text>
  <xsl:text>&lt;/td&gt;&#10;</xsl:text>
</xsl:template>

<xsl:template name="index-row">
  <xsl:param name="n"/>
  <xsl:for-each select="following-sibling::*[position() = 1]">
    <xsl:choose>
      <xsl:when test="not(@context) or @context=$hlp">
        <xsl:call-template name="index-cell">
          <xsl:with-param name="entry" select="@name"/>
        </xsl:call-template>        
        <xsl:if test="$n - 1 &gt; 0">
          <xsl:call-template name="index-row">
            <xsl:with-param name="n">
              <xsl:value-of select="$n - 1"/>
            </xsl:with-param>
          </xsl:call-template>
        </xsl:if>
      </xsl:when>
      <xsl:otherwise>
        <xsl:if test="$n &gt; 0">
          <xsl:call-template name="index-row">
            <xsl:with-param name="n">
              <xsl:value-of select="$n"/>
            </xsl:with-param>
          </xsl:call-template>
        </xsl:if>
      </xsl:otherwise>
    </xsl:choose>
  </xsl:for-each>
</xsl:template>

<xsl:template match="commandlist">
  <xsl:text>&lt;h1&gt;</xsl:text>
  <xsl:text>Gretl Command Reference</xsl:text>
  <xsl:text>&lt;/h1&gt;&#10;</xsl:text>
  <xsl:text>&lt;p&gt;</xsl:text>
  <xsl:text>See also the &lt;a href="./funcref.html"&gt;Gretl Function Reference&lt;/a&gt;.</xsl:text>
  <xsl:text>&lt;/p&gt;&#10;</xsl:text>
  <xsl:text>&lt;p&gt;</xsl:text>
  <xsl:text>The following commands are documented below.</xsl:text>
  <xsl:text>&lt;/p&gt;&#10;</xsl:text>
  <xsl:text>&lt;table cellspacing="4"&gt;&#10;</xsl:text>
  <xsl:for-each select="command[not(@context) or @context=$hlp]">
    <xsl:if test="position() mod 8 = 1">
      <xsl:text>&lt;tr&gt;&#10;</xsl:text>
      <xsl:call-template name="index-cell">
        <xsl:with-param name="entry" select="@name"/>
      </xsl:call-template>
      <xsl:call-template name="index-row">
        <xsl:with-param name="n" select="7"/>
      </xsl:call-template>
      <xsl:text>&lt;/tr&gt;&#10;</xsl:text>
    </xsl:if>
  </xsl:for-each>
  <xsl:text>&lt;/table&gt;&#10;</xsl:text>
  <xsl:text>&lt;p&gt;</xsl:text>
  <xsl:text>Note that brackets "&lt;kbd&gt;[&lt;/kbd&gt;" and "&lt;kbd&gt;]&lt;/kbd&gt;" </xsl:text>
  <xsl:text>are used to indicate that certain elements of commands are optional. </xsl:text>
  <xsl:text>The brackets should not be typed by the user.</xsl:text>
  <xsl:text>&lt;/p&gt;&#10;</xsl:text>
  <xsl:apply-templates/>
</xsl:template>

<xsl:template name="divider">
  <xsl:call-template name="nl"/>
  <xsl:text>&lt;hr&gt;</xsl:text>
  <xsl:call-template name="nl"/>
</xsl:template>

<xsl:template match="command">
  <xsl:if test="not(@context) or @context=$hlp">
    <xsl:if test="position() > 1">
      <xsl:call-template name="divider"/>
    </xsl:if>
    <xsl:text>&lt;h2&gt;</xsl:text>
    <xsl:text>&lt;a name="</xsl:text>
    <xsl:value-of select="@name"/>
    <xsl:text>"&gt;</xsl:text>
    <xsl:value-of select="@name"/>
    <xsl:text>&lt;/h2&gt;</xsl:text>
    <xsl:text>&#10;</xsl:text>
    <xsl:apply-templates/>
    <xsl:call-template name="dnl"/>
  </xsl:if>
</xsl:template>

<xsl:template match="description">
  <xsl:if test="not(@context) or @context=$hlp">
    <xsl:apply-templates/>
  </xsl:if>
</xsl:template>

<xsl:template match="usage">
  <xsl:if test="$hlp='cli'">
    <xsl:text>&lt;table cellspacing="4"&gt;</xsl:text>
    <xsl:apply-templates/>
    <xsl:text>&lt;/table&gt;</xsl:text>
    <xsl:call-template name="dnl"/>
  </xsl:if>
</xsl:template>

<xsl:template name="bold-cell-start">
  <xsl:text>&lt;td&gt;&lt;b&gt;</xsl:text>
</xsl:template>

<xsl:template name="bold-cell-end">
  <xsl:text>&lt;/b&gt;&lt;/td&gt;</xsl:text>
</xsl:template>

<xsl:template match="arguments">
  <xsl:text>&#10;&lt;tr&gt;</xsl:text>
  <xsl:call-template name="bold-cell-start"/>
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
  <xsl:call-template name="bold-cell-end"/>
  <xsl:text>&lt;td&gt;</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>&lt;/td&gt;</xsl:text>
  <xsl:text>&lt;/tr&gt;</xsl:text>
</xsl:template>

<xsl:template match="argblock">
  <xsl:if test="(@optional)">[ </xsl:if>
  <xsl:if test="(@separated)">; </xsl:if>
  <xsl:apply-templates/>
  <xsl:if test="(@optional)"> ] </xsl:if>
</xsl:template>

<xsl:template match="argument">
  <xsl:if test="(@optional)">[ </xsl:if> 
  <xsl:if test="(@separated)">; </xsl:if>
  <xsl:if test="(@alternate)">
    <xsl:call-template name="gettext">
      <xsl:with-param name="key" select="'or'"/>
    </xsl:call-template>
  </xsl:if>
  <xsl:if test="@flag">
    <xsl:text>&lt;kbd&gt;</xsl:text>
    <xsl:value-of select="@flag"/>
    <xsl:text>&lt;/kbd&gt;</xsl:text>
  </xsl:if> 
  <xsl:text>&lt;i&gt;</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>&lt;/i&gt;</xsl:text>
  <xsl:if test="(@optional)"> ] </xsl:if>
  <xsl:text> </xsl:text>
</xsl:template>

<xsl:template match="argpunct">
  <xsl:text>&lt;kbd&gt;</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>&lt;/kbd&gt;</xsl:text>
</xsl:template>

<xsl:template match="options">
  <xsl:call-template name="nl"/>
  <xsl:text>&lt;tr&gt;&lt;td&gt;&lt;b&gt;</xsl:text>
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
  <xsl:text>&lt;/b&gt;&lt;/td&gt;</xsl:text>
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="option">
  <xsl:if test="position() > 1">
    <xsl:text>&#10;&lt;tr&gt;&lt;td&gt;&lt;/td&gt;</xsl:text>
  </xsl:if>
  <xsl:text>&lt;td&gt;</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>&lt;/td&gt;&lt;/tr&gt;</xsl:text>
</xsl:template>

<xsl:template match="optparm">
  <xsl:if test="(@optional)">[</xsl:if>  
  <xsl:text>=&lt;i&gt;</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>&lt;/i&gt;</xsl:text>
  <xsl:if test="(@optional)">]</xsl:if>  
</xsl:template>

<xsl:template match="example">
  <xsl:if test="position() > 1">
    <xsl:text>&#10;&lt;tr&gt;&lt;td&gt;&lt;/td&gt;</xsl:text>
  </xsl:if> 
  <xsl:text>&lt;td&gt;&lt;kbd&gt;</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>&lt;/kbd&gt;&lt;/td&gt;&lt;/tr&gt;</xsl:text>
</xsl:template>

<xsl:template match="demos">
  <xsl:if test="preceding-sibling::*">
    <xsl:text>&lt;tr&gt;&lt;td&gt;&lt;/td&gt;</xsl:text>
  </xsl:if>
  <xsl:text>&lt;td&gt;</xsl:text>
  <xsl:if test="position() > 1">
    <xsl:call-template name="gettext">
      <xsl:with-param name="key" select="'Seealso'"/>
    </xsl:call-template>
  </xsl:if>
  <xsl:apply-templates/>
  <xsl:text>&lt;/td&gt;&lt;/tr&gt;</xsl:text>
  <xsl:call-template name="nl"/>
</xsl:template>

<xsl:template match="demo">
  <xsl:if test="position() > 1">
    <xsl:text>, </xsl:text>
  </xsl:if>
  <xsl:text>&lt;a href="scripts/</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>"&gt;</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>&lt;/a&gt;</xsl:text>
</xsl:template>

<xsl:template match="altforms">
  <xsl:call-template name="nl"/>
  <xsl:text>&lt;tr&gt;&lt;td&gt;&lt;b&gt;</xsl:text>
  <xsl:call-template name="gettext">
    <xsl:with-param name="key" select="'altforms'"/>
  </xsl:call-template>
  <xsl:text>&lt;/b&gt;&lt;/td&gt;</xsl:text>
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="altform">
  <xsl:if test="position() > 1">
    <xsl:text>&#10;&lt;tr&gt;&lt;td&gt;&lt;/td&gt;</xsl:text>
  </xsl:if> 
  <xsl:text>&lt;td&gt;</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>&lt;/td&gt;&lt;/tr&gt;</xsl:text>
</xsl:template>

<xsl:template match="syntax">
  <xsl:call-template name="nl"/>
  <xsl:text>&lt;tr&gt;&lt;td&gt;&lt;b&gt;</xsl:text>
  <xsl:call-template name="gettext">
    <xsl:with-param name="key" select="'syntax'"/>
  </xsl:call-template>
  <xsl:text>&lt;/b&gt;&lt;/td&gt;</xsl:text>
  <xsl:text>&lt;td&gt;</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>&lt;/td&gt;&lt;/tr&gt;</xsl:text>
</xsl:template>

<xsl:template match="examples">
  <xsl:if test="ancestor::command">
    <xsl:call-template name="nl"/>
  </xsl:if>
  <xsl:text>&lt;tr&gt;&lt;td&gt;&lt;b&gt;</xsl:text>
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
  <xsl:text>&lt;/b&gt;&lt;/td&gt;</xsl:text>
  <xsl:apply-templates/>
  <xsl:if test="ancestor::function">
    <xsl:call-template name="nl"/>
  </xsl:if>
</xsl:template>

<xsl:template match="effect">
  <xsl:text> (</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>)</xsl:text>
</xsl:template>

<xsl:template match="note">
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="funcref">
  <xsl:text>&lt;h1&gt;</xsl:text>
  <xsl:text>Gretl Function Reference</xsl:text>
  <xsl:text>&lt;/h1&gt;&#10;</xsl:text>
  <xsl:text>&lt;p&gt;</xsl:text>
  <xsl:text>See also the &lt;a href="./cmdref.html"&gt;Gretl Command Reference&lt;/a&gt;</xsl:text>
  <xsl:text>&#10;&lt;/p&gt;</xsl:text>
  <xsl:text>The following accessors and functions are documented below.</xsl:text>
  <xsl:text>&lt;/p&gt;&#10;</xsl:text>
  <xsl:for-each select="funclist">
    <xsl:text>&lt;table cellspacing="4"&gt;&#10;</xsl:text>
    <xsl:for-each select="function">
      <xsl:if test="position() mod 8 = 1">
        <xsl:text>&lt;tr&gt;&#10;</xsl:text>
        <xsl:call-template name="index-cell">
          <xsl:with-param name="entry" select="@name"/>
        </xsl:call-template>
        <xsl:call-template name="index-row">
          <xsl:with-param name="n" select="7"/>
        </xsl:call-template>
        <xsl:text>&lt;/tr&gt;&#10;</xsl:text>
      </xsl:if>
    </xsl:for-each>
    <xsl:text>&lt;/table&gt;&#10;</xsl:text>
  </xsl:for-each>
  <xsl:apply-templates/> 
</xsl:template>

<xsl:template match="funclist">
  <xsl:text>&#10;&lt;h1&gt;</xsl:text>
  <xsl:value-of select="@name"/>
  <xsl:text>&lt;/h1&gt;&#10;&#10;</xsl:text>
  <xsl:apply-templates/> 
</xsl:template>

<xsl:template match="function">
  <xsl:if test="position() > 1">
    <xsl:call-template name="divider"/>
  </xsl:if>
  <xsl:text>&lt;h2&gt;</xsl:text>
  <xsl:text>&lt;a name="</xsl:text>
  <xsl:value-of select="@name"/>
  <xsl:text>"&gt;</xsl:text>
  <xsl:value-of select="@name"/>
  <xsl:text>&lt;/h2&gt;&#10;</xsl:text>
  <xsl:text>&lt;b&gt;</xsl:text>
  <xsl:call-template name="gettext">
    <xsl:with-param name="key" select="'output'"/>
  </xsl:call-template>
  <xsl:text>&lt;/b&gt; </xsl:text>
  <xsl:call-template name="gettext">
    <xsl:with-param name="key" select="@output"/>
  </xsl:call-template>
  <xsl:if test="@altout">
    <xsl:call-template name="gettext">
      <xsl:with-param name="key" select="'or'"/>
    </xsl:call-template>
    <xsl:call-template name="gettext">
      <xsl:with-param name="key" select="@altout"/>
    </xsl:call-template>
  </xsl:if>
  <xsl:if test="not(@fnargs)">
    <xsl:text>&#10;</xsl:text>
  </xsl:if>
  <xsl:apply-templates/>
  <xsl:call-template name="dnl"/>
</xsl:template>

<xsl:template match="fnargs">
  <xsl:text>&lt;table&gt;</xsl:text>
  <xsl:text>&lt;tr&gt;&lt;td&gt;&lt;b&gt;</xsl:text>
  <xsl:choose>
    <xsl:when test="count(fnarg) > 1">
      <xsl:call-template name="gettext-nospace">
        <xsl:with-param name="key" select="'args'"/>
      </xsl:call-template>
    </xsl:when>
    <xsl:otherwise>
      <xsl:call-template name="gettext">
        <xsl:with-param name="key" select="'arg'"/>
      </xsl:call-template>
    </xsl:otherwise> 
  </xsl:choose>
  <xsl:text>&lt;/b&gt;&lt;/td&gt;</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>&lt;/table&gt;</xsl:text>
</xsl:template>

<xsl:template match="fnarg">
  <xsl:if test="position() > 1">
    <xsl:text>&lt;tr&gt;&lt;td&gt;&lt;/td&gt;</xsl:text>
  </xsl:if>
  <xsl:text>&lt;td&gt;</xsl:text>
  <xsl:choose>
    <xsl:when test="@type='varargs'">
      <xsl:text>...</xsl:text>
    </xsl:when>
    <xsl:otherwise>
      <xsl:text>&lt;i&gt;</xsl:text>
      <xsl:apply-templates/>
      <xsl:text>&lt;/i&gt;</xsl:text>
    </xsl:otherwise>
  </xsl:choose>
  <xsl:text> (</xsl:text>
  <xsl:call-template name="gettext-nospace">
    <xsl:with-param name="key" select="@type"/>
  </xsl:call-template>
  <xsl:if test="(@optional)">
    <xsl:choose>
      <xsl:when test="@type='matrixref'">
        <xsl:text>, </xsl:text>
        <xsl:call-template name="gettext">
          <xsl:with-param name="key" select="'or'"/>
        </xsl:call-template>
        <xsl:text>&lt;kbd&gt;null&lt;/kbd&gt;</xsl:text>
      </xsl:when>
      <xsl:otherwise>
        <xsl:text>, </xsl:text>
        <xsl:call-template name="gettext-nospace">
          <xsl:with-param name="key" select="'optional'"/>
        </xsl:call-template>
      </xsl:otherwise>
    </xsl:choose>
  </xsl:if> 
  <xsl:text>)</xsl:text>
  <xsl:text>&lt;/td&gt;&lt;/tr&gt;&#10;</xsl:text>
</xsl:template>

<xsl:template match="repl|argname">
  <xsl:text>&lt;i&gt;</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>&lt;/i&gt;</xsl:text>
</xsl:template>

<xsl:template match="emphasis">
  <xsl:text>&lt;i&gt;</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>&lt;/i&gt;</xsl:text>
</xsl:template>

<xsl:template match="cmd">
  <xsl:text>&lt;kbd&gt;</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>&lt;/kbd&gt;</xsl:text>
</xsl:template>

<xsl:template match="opt">
  <xsl:text>&lt;kbd&gt;--</xsl:text>
  <xsl:choose>
    <xsl:when test="substring(text(),1,2)='--'">
      <xsl:value-of select='substring-after(text(),"--")'/>
    </xsl:when>
    <xsl:otherwise>
      <xsl:apply-templates/>
    </xsl:otherwise>
  </xsl:choose>
  <xsl:text>&lt;/kbd&gt;</xsl:text>
</xsl:template>

<xsl:template match="program">
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="lit|func|filename">
  <xsl:text>&lt;kbd&gt;</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>&lt;/kbd&gt;</xsl:text>
</xsl:template>

<xsl:template match="flag">
  <xsl:text>&lt;kbd&gt;</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>&lt;/kbd&gt;</xsl:text>
</xsl:template>

<xsl:template match="math">
  <xsl:text>&lt;i&gt;</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>&lt;/i&gt;</xsl:text>
</xsl:template>

<xsl:template match="sup">
  <xsl:text>&lt;sup&gt;</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>&lt;/sup&gt;</xsl:text>
</xsl:template>

<xsl:template match="sub">
  <xsl:text>&lt;sub&gt;</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>&lt;/sub&gt;</xsl:text>
</xsl:template>

<xsl:template match="book">
  <xsl:text>&lt;i&gt;</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>&lt;/i&gt;</xsl:text>
</xsl:template>

<xsl:template match="cite">
  <xsl:text>&lt;a href="./biblio.html#</xsl:text>
  <xsl:value-of select="@key"/>
  <xsl:text>"&gt;</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>&lt;/a&gt;</xsl:text>
</xsl:template>

<xsl:template match="quote">
  <xsl:text>"</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>"</xsl:text>
</xsl:template>

<xsl:template match="url">
  <xsl:text>&lt;a href="http://</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>"&gt;</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>&lt;/a&gt;</xsl:text>
</xsl:template>

<xsl:template match="equation">
  <xsl:if test="(@status='display')">
    <xsl:text>&lt;/p&gt;&#10;&#10;</xsl:text>
  </xsl:if>
  <xsl:choose>
    <xsl:when test="@graphic">
      <xsl:text>&lt;img alt="missing image" src="./images/</xsl:text>
      <xsl:value-of select="@graphic"/>
      <xsl:text>.png"&gt;</xsl:text>
    </xsl:when>
    <xsl:otherwise>
      <xsl:text>&lt;i&gt;</xsl:text>
      <xsl:value-of select="@ascii"/>
      <xsl:text>&lt;/i&gt;</xsl:text>
    </xsl:otherwise>
  </xsl:choose>
  <xsl:if test="(@status='display')">
    <xsl:text>&#10;&#10;&lt;p&gt;</xsl:text>
  </xsl:if>
</xsl:template>

<xsl:template match="para">
  <xsl:if test="not(@context) or @context=$hlp or @context='notex'">
    <xsl:text>&#10;&lt;p&gt;</xsl:text>
    <xsl:apply-templates/>
    <xsl:text>&lt;/p&gt;</xsl:text>
 </xsl:if>
</xsl:template>

<xsl:template match="gfr">
  <xsl:text>the &lt;a href="./funcref.html"&gt;</xsl:text>
  <xsl:text>Gretl Function Reference</xsl:text>
  <xsl:text>&lt;/a&gt;</xsl:text>
</xsl:template>

<xsl:template match="refnote">
  <xsl:if test="@xref='true'">
    <xsl:apply-templates/>
  </xsl:if>
</xsl:template>

<xsl:template match="code">
  <xsl:if test="not(@context) or @context=$hlp or @context='notex'">
    <xsl:text>&#10;&lt;pre&gt;</xsl:text>
    <xsl:apply-templates/>
    <xsl:text>&lt;/pre&gt;&#10;</xsl:text>
  </xsl:if>
</xsl:template>

<xsl:template match="mono">
  <xsl:if test="not(@context) or @context=$hlp or @context='notex'">
    <xsl:text>&lt;pre&gt;</xsl:text>
    <xsl:apply-templates/>
    <xsl:text>&lt;/pre&gt;</xsl:text>
  </xsl:if>
</xsl:template>

<xsl:template match="pre">
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="ilist">
  <xsl:if test="not(@context) or @context=$hlp or @context='notex'">
    <xsl:text>&lt;ul&gt;</xsl:text>
    <xsl:apply-templates/>
    <xsl:text>&lt;/ul&gt;</xsl:text>
  </xsl:if>
</xsl:template>

<xsl:template match="nlist">
  <xsl:if test="not(@context) or @context=$hlp">
    <xsl:text>&lt;ol&gt;</xsl:text>
    <xsl:apply-templates/>
    <xsl:text>&lt;/ol&gt;</xsl:text>
  </xsl:if>
</xsl:template>

<xsl:template match="li">
  <xsl:text>&lt;li&gt;</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>&lt;/li&gt;</xsl:text>
</xsl:template>

<xsl:template match="subhead">
  <xsl:if test="not(@context) or @context=$hlp">
    <xsl:text>&lt;h3&gt;</xsl:text>
    <xsl:apply-templates/>
    <xsl:text>&lt;/h3&gt;</xsl:text>
  </xsl:if>
</xsl:template>

<xsl:template match="cmdref">
  <xsl:choose>
    <xsl:when test="$topic = 'funcs'">
      <xsl:text>&lt;a href="./cmdref.html#</xsl:text>
    </xsl:when>
    <xsl:otherwise>
      <xsl:text>&lt;a href="#</xsl:text>
    </xsl:otherwise>
  </xsl:choose>
  <xsl:value-of select="@targ"/>
  <xsl:text>"&gt;</xsl:text>
  <xsl:value-of select="@targ"/>
  <xsl:text>&lt;/a&gt;</xsl:text>
  <xsl:if test="parent::seelist and following-sibling::*">
    <xsl:text>, </xsl:text>
  </xsl:if>
</xsl:template>

<xsl:template match="fncref">
  <xsl:choose>
    <xsl:when test="$topic = 'funcs'">
      <xsl:text>&lt;a href="#</xsl:text>
    </xsl:when>
    <xsl:otherwise>
      <xsl:text>&lt;a href="./funcref.html#</xsl:text>
    </xsl:otherwise>
  </xsl:choose>
  <xsl:choose>
    <xsl:when test="@label">
      <xsl:value-of select="substring(@label, 2)"/>
    </xsl:when>
    <xsl:otherwise>
      <xsl:value-of select="@targ"/>
    </xsl:otherwise>
  </xsl:choose>
  <xsl:text>"&gt;</xsl:text>
  <xsl:choose>
    <xsl:when test="@label">
      <xsl:value-of select="substring(@label, 2)"/>
    </xsl:when>
    <xsl:otherwise>
      <xsl:value-of select="@targ"/>
    </xsl:otherwise>
  </xsl:choose>
  <xsl:text>&lt;/a&gt;</xsl:text>
  <xsl:if test="parent::seelist and following-sibling::*">
    <xsl:text>, </xsl:text>
  </xsl:if>
</xsl:template>

<xsl:template match="guideref">
  <xsl:text>&lt;a href="gretl-guide.pdf"&gt;</xsl:text>
  <xsl:call-template name="gettext-nospace">
    <xsl:with-param name="key" select="'guidebook'"/>
  </xsl:call-template>
  <xsl:text>&lt;/a&gt;</xsl:text>
</xsl:template>

<xsl:template match="menu-path">
  <xsl:text>&#10;&lt;p&gt;</xsl:text>
  <xsl:text>&lt;b&gt;</xsl:text>
  <xsl:call-template name="gettext">
    <xsl:with-param name="key" select="'menupath'"/>
  </xsl:call-template>
  <xsl:text>&lt;/b&gt;</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>&lt;/p&gt;</xsl:text>
</xsl:template>

<xsl:template match="other-access">
  <xsl:text>&#10;</xsl:text>
  <xsl:call-template name="gettext">
    <xsl:with-param name="key" select="'otheraccess'"/>
  </xsl:call-template>
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="table"> 
<!-- can't handle tables at present -->
</xsl:template>

<xsl:template match="footnote"/>

<xsl:template name="nl">
  <xsl:text>&#10;</xsl:text>  
</xsl:template>

<xsl:template name="dnl">
  <xsl:text>&#10;&#10;</xsl:text>  
</xsl:template>

<xsl:template match="seelist">
  <xsl:call-template name="gettext">
    <xsl:with-param name="key" select="'Seealso'"/>
  </xsl:call-template>
  <xsl:apply-templates/>
  <xsl:text>.</xsl:text>
</xsl:template>

<xsl:template match="by">
  <xsl:choose>
    <xsl:when test="contains(@r,'1') or contains(@r,'2')">
      <xsl:value-of select="@r"/>
    </xsl:when>    
    <xsl:otherwise>
      <xsl:text>&lt;i&gt;</xsl:text> 
      <xsl:value-of select="@r"/> 
      <xsl:text>&lt;/i&gt;</xsl:text>
    </xsl:otherwise>
  </xsl:choose>
  <xsl:text> x </xsl:text> 
  <xsl:choose>
    <xsl:when test="contains(@r,'1') or contains(@r,'2')">
      <xsl:value-of select="@c"/>
    </xsl:when>    
    <xsl:otherwise>
      <xsl:text>&lt;i&gt;</xsl:text> 
      <xsl:value-of select="@c"/> 
      <xsl:text>&lt;/i&gt;</xsl:text> 
    </xsl:otherwise>
  </xsl:choose>
</xsl:template>

</xsl:stylesheet>
