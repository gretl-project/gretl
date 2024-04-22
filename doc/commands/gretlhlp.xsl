<?xml version="1.0"?>
<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform" >

  <!-- Stylesheet for formatted GUI "online" help -->

<xsl:param name="hlp">cli</xsl:param>
<xsl:param name="fmt">plain</xsl:param>  
<xsl:param name="topic">cmds</xsl:param>
<xsl:param name="refs">chaprefs.xml</xsl:param>
<xsl:param name="lang" select="'en'"/>

<xsl:output method="text" encoding="utf-8"/>

<xsl:variable name="intl"
	      select="document('hlp_l10n.xml')/internationalization"/>
  
<xsl:variable name="docref"
	      select="document($refs)/refsets"/>

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

<xsl:template name="getref">
  <xsl:param name="key"/>
  <xsl:variable name="chap"
    select="normalize-space($docref/refset[@id='guide-chapters']/ref[@key=$key]/@chapter)"/>
  <xsl:choose>
    <xsl:when test="$chap">
      <xsl:value-of select="$chap"/>
    </xsl:when>
    <xsl:otherwise>
      <xsl:message terminate="yes">
        <xsl:text>** Error: no phrase with key = '</xsl:text>
        <xsl:value-of select="$key"/>
        <xsl:text>' found.'</xsl:text>
      </xsl:message>
    </xsl:otherwise>
  </xsl:choose>  
</xsl:template>

<xsl:template match="commandref">
  <xsl:apply-templates/> 
</xsl:template>

<xsl:template match="command">
  <xsl:if test="not(@context) or @context=$hlp">
    <xsl:if test="position() > 1">
      <xsl:call-template name="nl"/>
    </xsl:if>
    <xsl:text># </xsl:text>
    <xsl:value-of select="@name"/>
    <xsl:text> </xsl:text>
    <xsl:value-of select="@section"/>
    <xsl:if test="$hlp='gui'">
      <xsl:text> "</xsl:text>
      <xsl:value-of select="@label"/>
      <xsl:text>"</xsl:text>
    </xsl:if>
    <xsl:text>&#10;</xsl:text>
    <xsl:apply-templates/>
    <xsl:call-template name="dnl"/>
    <xsl:if test="(not(@context) and $hlp='gui')">
      <xsl:call-template name="gettext">
	<xsl:with-param name="key" select="'scriptcommand'"/>
      </xsl:call-template>
      <xsl:text>&lt;@ref="</xsl:text>
      <xsl:value-of select="@name"/>
      <xsl:text>"&gt;&#10;</xsl:text>
    </xsl:if>
  </xsl:if>
</xsl:template>

<xsl:template match="common-opt" />

<xsl:template match="not-ready-common-opt">
  <xsl:if test="position() > 1">
    <xsl:call-template name="nl"/>
  </xsl:if>
  <xsl:text># option:</xsl:text>
  <xsl:value-of select="@name"/>
  <xsl:call-template name="dnl"/>
  <xsl:text>&lt;@lit="--</xsl:text>
  <xsl:value-of select="@name"/>
  <xsl:if test="@param">
    <xsl:text>=</xsl:text>
    <xsl:text>"&gt;</xsl:text>
    <xsl:text>&lt;@var="</xsl:text>
    <xsl:value-of select="@param"/>
  </xsl:if>
  <xsl:text>"&gt;</xsl:text>
  <xsl:text>&#10;</xsl:text>
  <xsl:apply-templates/>
  <xsl:call-template name="dnl"/>
</xsl:template>

<xsl:template match="applies-to">
  <xsl:if test="position() = 1">
    <xsl:text>&#10;</xsl:text>
    <xsl:call-template name="gettext-nospace">
      <xsl:with-param name="key" select="'applies-to'"/>
    </xsl:call-template>
    <xsl:text>:&#10;</xsl:text>
  </xsl:if>
  <xsl:if test="position() > 1">
    <xsl:call-template name="nl"/>
  </xsl:if>
  <xsl:text>&lt;@lit="</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>"&gt;</xsl:text>
  <xsl:if test="not(following-sibling::applies-to)">
    <xsl:text>&#10;</xsl:text>
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
  <xsl:text>&#x9;</xsl:text>
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="argblock">
  <xsl:if test="(@optional)">[ </xsl:if>
  <xsl:if test="(@separated)">; </xsl:if>
  <xsl:apply-templates/>
  <xsl:if test="(@optional)">] </xsl:if>
</xsl:template>

<xsl:template match="argument">
  <xsl:if test="(@optional)">[ </xsl:if> 
  <xsl:if test="(@separated and (preceding-sibling::argument or preceding-sibling::argblock))">; </xsl:if>
  <xsl:if test="(@alternate)">
    <xsl:call-template name="gettext">
      <xsl:with-param name="key" select="'or'"/>
    </xsl:call-template>
  </xsl:if>
  <xsl:if test="@flag">
    <xsl:text>&lt;@lit="</xsl:text>
    <xsl:value-of select="@flag"/>
    <xsl:text>"&gt;</xsl:text>
  </xsl:if> 
  <xsl:text>&lt;@var="</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>"&gt; </xsl:text>
  <xsl:if test="(@separated and not(preceding-sibling::argument or preceding-sibling::argblock))">; </xsl:if>
  <xsl:if test="(@optional)">] </xsl:if> 
</xsl:template>

<xsl:template match="argpunct">
  <xsl:text>&lt;@lit="</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>"&gt;</xsl:text>
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

<xsl:template match="option">
  <xsl:if test="position() > 1">
    <xsl:text>&#xa;&#x9;</xsl:text>
  </xsl:if>
  <xsl:text>&#x9;</xsl:text>
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="optnote">
  <xsl:text>&#xa;&#x9;</xsl:text>
  <xsl:text>&#x9;</xsl:text>
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="optparm">
  <xsl:if test="(@optional)">[</xsl:if>  
  <xsl:text>=&lt;@var="</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>"&gt;</xsl:text>
  <xsl:if test="(@optional)">]</xsl:if>  
</xsl:template>

<xsl:template match="example">
  <xsl:if test="position() > 1">
    <xsl:text>&#xa;&#x9;</xsl:text>
  </xsl:if> 
  <xsl:text>&#x9;</xsl:text>
  <xsl:text>&lt;@lit="</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>"&gt;</xsl:text>  
</xsl:template>

<xsl:template match="demos">
  <xsl:choose>
    <xsl:when test="position() > 1">
      <xsl:text>&#xa;&#x9;&#x9;</xsl:text>
      <xsl:call-template name="gettext">
        <xsl:with-param name="key" select="'Seealso'"/>
      </xsl:call-template>
    </xsl:when>
    <xsl:otherwise>
      <xsl:text>&#x9;</xsl:text>
    </xsl:otherwise>
  </xsl:choose>
  <xsl:apply-templates/>
  <xsl:call-template name="nl"/>
</xsl:template>

<xsl:template match="demo">
  <xsl:if test="position() > 1">
    <xsl:text>, </xsl:text>
  </xsl:if>
  <xsl:text>&lt;@inp="</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>"&gt;</xsl:text>
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
    <xsl:text>&#xa;&#x9;</xsl:text>
  </xsl:if> 
  <xsl:text>&#x9;</xsl:text>
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="syntax">
  <xsl:call-template name="nl"/>
  <xsl:call-template name="gettext">
    <xsl:with-param name="key" select="'syntax'"/>
  </xsl:call-template>
  <xsl:text>&#x9;</xsl:text>
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="examples">
  <xsl:if test="ancestor::command">
    <xsl:call-template name="nl"/>
  </xsl:if>
  <xsl:choose>
    <xsl:when test="count(example) > 1 or count(demos) > 0">
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
  <xsl:if test="ancestor::function">
    <xsl:call-template name="nl"/>
  </xsl:if>
</xsl:template>

<xsl:template match="effect">
  <xsl:text> (</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>)</xsl:text>
</xsl:template>

<xsl:template match="funcref">
  <xsl:apply-templates/> 
</xsl:template>

<xsl:template match="funclist">
  <xsl:text>&#10;</xsl:text>
  <xsl:text>## </xsl:text>
  <xsl:value-of select="@name"/>
  <xsl:text>&#10;&#10;</xsl:text>
  <xsl:apply-templates/> 
</xsl:template>

<xsl:template match="function">
  <xsl:if test="position() > 1">
    <xsl:call-template name="nl"/>
  </xsl:if>
  <xsl:text># </xsl:text>
  <xsl:value-of select="@name"/>
  <xsl:text> </xsl:text>
  <xsl:value-of select="@section"/>
  <xsl:text>&#10;</xsl:text>
  <xsl:call-template name="gettext">
    <xsl:with-param name="key" select="'output'"/>
  </xsl:call-template>
  <xsl:text>&#x9;</xsl:text>
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
  <xsl:text>&#x9;</xsl:text>
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="fnarg">
  <xsl:if test="position() > 1">
    <xsl:text>&#x9;&#x9;</xsl:text>
  </xsl:if> 
  <xsl:choose>
    <xsl:when test="@type='varargs'">
      <xsl:text>...</xsl:text>
    </xsl:when>
    <xsl:otherwise>
      <xsl:text>&lt;@var="</xsl:text>
      <xsl:apply-templates/>
      <xsl:text>"&gt; </xsl:text>
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
        <xsl:text>&lt;@lit="null"&gt;</xsl:text>
      </xsl:when>
      <xsl:otherwise>
        <xsl:text>, </xsl:text>
        <xsl:call-template name="gettext-nospace">
          <xsl:with-param name="key" select="'optional'"/>
        </xsl:call-template>
      </xsl:otherwise>
    </xsl:choose>
  </xsl:if>
  <xsl:if test="(@conditional)">
    <xsl:text>, </xsl:text>
    <xsl:call-template name="gettext-nospace">
      <xsl:with-param name="key" select="'conditional'"/>
    </xsl:call-template>
  </xsl:if> 
  <xsl:text>)</xsl:text>
  <xsl:text>&#10;</xsl:text>
</xsl:template>

<xsl:template match="repl|argname">
  <xsl:text>&lt;@var="</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>"&gt;</xsl:text>
</xsl:template>

<xsl:template match="emphasis|i">
  <xsl:text>&lt;@itl="</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>"&gt;</xsl:text>
</xsl:template>

<xsl:template match="cmd">
  <xsl:text>&lt;@lit="</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>"&gt;</xsl:text>
</xsl:template>

<xsl:template match="opt">
  <xsl:comment>
    For option flags, ensure that there's no break between 
    the double-dash and the following string.
  </xsl:comment>
  <xsl:text>&lt;@opt="--&#x2060;</xsl:text>
  <xsl:choose>
    <xsl:when test="substring(text(),1,2)='--'">
      <xsl:value-of select='substring-after(text(),"--")'/>
    </xsl:when>
    <xsl:otherwise>
      <xsl:apply-templates/>
    </xsl:otherwise>
  </xsl:choose>
  <xsl:text>"&gt;</xsl:text>
</xsl:template>

<xsl:template match="program">
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="lit|func|filename|flag">
  <xsl:text>&lt;@lit="</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>"&gt;</xsl:text>
</xsl:template>

<xsl:template match="math">
  <xsl:text>&lt;@mth="</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>"&gt;</xsl:text>
</xsl:template>

<xsl:template match="sup">
  <xsl:text>&lt;@sup="</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>"&gt;</xsl:text>
</xsl:template>

<xsl:template match="sub">
  <xsl:text>&lt;@sub="</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>"&gt;</xsl:text>
</xsl:template>

<xsl:template match="book">
  <xsl:text>&lt;@itl="</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>"&gt;</xsl:text>
</xsl:template>

<xsl:template match="cite">
  <xsl:text>&lt;@bib="</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>;</xsl:text>
  <xsl:value-of select="@key"/>
  <xsl:text>"&gt;</xsl:text>
</xsl:template>

<xsl:template match="quote">
  <xsl:text>&#x201c;</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>&#x201d;</xsl:text>
</xsl:template>

<xsl:template match="url">
  <xsl:choose>
    <xsl:when test="$fmt='pango'">
      <xsl:text>&lt;@url="</xsl:text>
    </xsl:when>
    <xsl:otherwise>
      <xsl:text>&lt;@lit="</xsl:text>
    </xsl:otherwise>
  </xsl:choose>
  <xsl:apply-templates/>
  <xsl:text>"&gt;</xsl:text>
</xsl:template>

<xsl:template match="doc">
  <xsl:choose>
    <xsl:when test="$fmt='pango'">
      <xsl:text>&lt;@adb="</xsl:text>
    </xsl:when>
    <xsl:otherwise>
      <xsl:text>&lt;@lit="</xsl:text>
    </xsl:otherwise>
  </xsl:choose>
  <xsl:apply-templates/>
  <xsl:text>"&gt;</xsl:text>
</xsl:template>

<xsl:template match="mnu">
  <xsl:choose>
    <xsl:when test="$fmt='pango'">
      <xsl:text>&lt;@mnu="</xsl:text>
      <xsl:value-of select="@targ"/>
      <xsl:text>"&gt;</xsl:text>
    </xsl:when>
    <xsl:otherwise>
      <xsl:apply-templates/>
    </xsl:otherwise>
  </xsl:choose>
</xsl:template>

<xsl:template match="equation">
  <xsl:if test="(@status='display')">
    <xsl:text>[/PARA]&#xa;&#xa;&#xa;  </xsl:text>
  </xsl:if>
  <xsl:choose>
    <xsl:when test="@graphic">
      <xsl:text>&lt;@fig="</xsl:text>
      <xsl:value-of select="@graphic"/>
      <xsl:text>"&gt;</xsl:text>
    </xsl:when>
    <xsl:otherwise>
      <xsl:text>&lt;@mth="</xsl:text>
      <xsl:value-of select="@ascii"/>
      <xsl:text>"&gt;</xsl:text>
    </xsl:otherwise>
  </xsl:choose>
  <xsl:if test="(@status='display')">
    <xsl:text>&#xa;&#xa;&#xa;[PARA]</xsl:text>
  </xsl:if>
</xsl:template>

<xsl:template match="para">
  <xsl:if test="not(@context) or @context=$hlp or @context='notex'">
  <xsl:choose>
    <xsl:when test="parent::li and ancestor::ilist">
      <xsl:text>&#xa;[ILISTPAR]</xsl:text>
      <xsl:if test="count(../*)=1">
        <xsl:text>&#x2022; </xsl:text>
      </xsl:if>
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

<xsl:template match="gfr">
  <xsl:text>&lt;@gfr="</xsl:text>
  <xsl:call-template name="gettext-nospace">
    <xsl:with-param name="key" select="'funcchap'"/>
  </xsl:call-template>  
  <xsl:text>"&gt;</xsl:text>
</xsl:template>

<xsl:template match="refnote">
  <xsl:if test="@xref='true'">
    <xsl:apply-templates/>
  </xsl:if>
</xsl:template>

<xsl:template match="code">
  <xsl:if test="not(@context) or @context=$hlp or @context='notex'">
    <xsl:text>[CODE]</xsl:text>
    <xsl:apply-templates/>
    <xsl:text>[/CODE]&#xa;</xsl:text>
  </xsl:if>
</xsl:template>

<xsl:template match="mono">
  <xsl:if test="not(@context) or @context=$hlp or @context='notex'">
    <xsl:text>[MONO]</xsl:text>
    <xsl:apply-templates/>
    <xsl:text>[/MONO]&#xa;</xsl:text>
  </xsl:if>
</xsl:template>

<xsl:template match="pre">
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="ilist">
  <xsl:if test="not(@context) or @context=$hlp or @context='notex'">
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

<xsl:template match="subhead">
  <xsl:if test="not(@context) or @context=$hlp">
    <xsl:text>[PARA]</xsl:text>
    <xsl:text>&lt;subhead&gt;</xsl:text>
    <xsl:apply-templates/>
    <xsl:text>&lt;/subhead&gt;</xsl:text>
    <xsl:text>[/PARA]</xsl:text>
  </xsl:if>
</xsl:template>

<xsl:template match="cmdref">
  <xsl:choose>
    <xsl:when test="$topic = 'funcs'">
      <xsl:text>&lt;@xrf="</xsl:text>
    </xsl:when>
    <xsl:otherwise>
      <xsl:text>&lt;@ref="</xsl:text>
    </xsl:otherwise>
  </xsl:choose>
  <xsl:value-of select="@targ"/>
  <xsl:text>"&gt;</xsl:text>
  <xsl:if test="parent::seelist and following-sibling::*">
    <xsl:text>, </xsl:text>
  </xsl:if>
</xsl:template>

<xsl:template match="fncref">
  <xsl:choose>
    <xsl:when test="$topic = 'funcs'">
      <xsl:text>&lt;@ref="</xsl:text>
    </xsl:when>
    <xsl:otherwise>
      <xsl:text>&lt;@xrf="</xsl:text>
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
  <xsl:if test="parent::seelist and following-sibling::*">
    <xsl:text>, </xsl:text>
  </xsl:if>
</xsl:template>

<xsl:template match="guideref">
  <xsl:choose>
    <xsl:when test="$lang = 'en'">
      <xsl:text>chapter </xsl:text>
      <xsl:call-template name="getref">
	<xsl:with-param name="key" select="@targ"/>
      </xsl:call-template>
      <xsl:text> of the &lt;@pdf="Gretl User's Guide#</xsl:text>
      <xsl:value-of select="@targ"/>
      <xsl:text>"&gt;</xsl:text>
    </xsl:when>
    <xsl:otherwise>
      <xsl:text>&lt;@pdf="</xsl:text>
      <xsl:call-template name="gettext-nospace">
	<xsl:with-param name="key" select="'guidebook'"/>
      </xsl:call-template>
      <xsl:text>#</xsl:text>
      <xsl:value-of select="@targ"/>
      <xsl:text>"&gt; (</xsl:text>
      <xsl:call-template name="gettext">
	<xsl:with-param name="key" select="'chapter'"/>
      </xsl:call-template>
      <xsl:call-template name="getref">
	<xsl:with-param name="key" select="@targ"/>
      </xsl:call-template>
      <xsl:text>)</xsl:text>
    </xsl:otherwise>
  </xsl:choose>
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
<!-- can't handle tables at present -->
</xsl:template>

<xsl:template match="tabular">
<!-- reserved for TeX -->
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
    <xsl:when test="number(@r) = @r">
      <xsl:value-of select="@r"/>
    </xsl:when>    
    <xsl:otherwise>
      <xsl:text>&lt;@itl="</xsl:text> 
      <xsl:value-of select="@r"/> 
      <xsl:text>"&gt;</xsl:text> 
    </xsl:otherwise>
  </xsl:choose>
  <xsl:text>&#215;</xsl:text> 
  <xsl:choose>
    <xsl:when test="number(@c) = @c">
      <xsl:value-of select="@c"/>
    </xsl:when>    
    <xsl:otherwise>
      <xsl:text>&lt;@itl="</xsl:text> 
      <xsl:value-of select="@c"/> 
      <xsl:text>"&gt;</xsl:text> 
    </xsl:otherwise>
  </xsl:choose>
</xsl:template>

</xsl:stylesheet>
