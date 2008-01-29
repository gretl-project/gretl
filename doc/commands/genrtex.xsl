<?xml version="1.0"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" >

<xsl:param name="hlp">tex</xsl:param>
<xsl:param name="xrefs">true</xsl:param>
<xsl:param name="useGUG">true</xsl:param>
<xsl:param name="lang" select="'en'"/>

<xsl:output method="text" encoding="iso-8859-1"/>

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

<xsl:strip-space elements="title"/>
<xsl:preserve-space elements="*"/>

<xsl:template match="funcref">
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="funclist|commandlist">
  <xsl:text>\section{</xsl:text>
  <xsl:value-of select="@name"/>
  <xsl:text>}&#10;\label{sec:</xsl:text>
  <xsl:value-of select="@ref"/>
  <xsl:text>}&#10;&#10;</xsl:text>
  <xsl:apply-templates/>
</xsl:template>

<xsl:template name="funcname">
  <xsl:param name="src"/>
  <xsl:if test="substring($src, 1, 1) = '$'">
    <xsl:text>\</xsl:text>
  </xsl:if>
  <xsl:value-of select="$src"/>  
</xsl:template>

<xsl:template name="escdollar">
  <xsl:param name="src"/>
  <xsl:choose>
    <xsl:when test="substring($src, 1, 1) = '$'">
      <xsl:text>\</xsl:text>
      <xsl:value-of select="$src"/>
    </xsl:when>
    <xsl:otherwise>
      <xsl:value-of select="$src"/>
    </xsl:otherwise>
  </xsl:choose>
</xsl:template>

<xsl:template match="function">
  <xsl:text>&#10;</xsl:text>
  <xsl:text>\subsection{</xsl:text>
  <xsl:call-template name="funcname">
    <xsl:with-param name="src" select="@name"/>
  </xsl:call-template>
  <xsl:text>}&#10;</xsl:text>
  <xsl:text>\hypertarget{func-</xsl:text>
  <xsl:choose>
    <xsl:when test="@targ">
      <xsl:value-of select="@targ"/>
    </xsl:when>
    <xsl:otherwise>
      <xsl:call-template name="escdollar">
        <xsl:with-param name="src" select="@name"/>
      </xsl:call-template>
    </xsl:otherwise>
  </xsl:choose>
  <xsl:text>}{}</xsl:text>
  <xsl:text>&#10;&#10;</xsl:text>
  <xsl:text>\begin{tabular}{ll}&#10;</xsl:text>
  <xsl:call-template name="gettext">
    <xsl:with-param name="key" select="'output'"/>
  </xsl:call-template>
  <xsl:text>&amp; </xsl:text>
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
  <xsl:text>\\&#10;</xsl:text>
  <xsl:if test="not(fnargs)">
    <xsl:text>\end{tabular}&#10;</xsl:text>
  </xsl:if>
  <xsl:apply-templates/>
  <xsl:call-template name="dnl"/>
</xsl:template>

<xsl:template match="command">
  <xsl:if test="@context != 'gui'">
    <xsl:text>&#10;</xsl:text>
    <xsl:text>\subsection{</xsl:text>
    <xsl:call-template name="funcname">
      <xsl:with-param name="src" select="@name"/>
    </xsl:call-template>
    <xsl:text>}&#10;</xsl:text>
    <xsl:text>\hypertarget{cmd-</xsl:text>
    <xsl:value-of select="@name"/>
    <xsl:text>}{}</xsl:text>
    <xsl:text>&#10;</xsl:text>
    <xsl:apply-templates/>
    <xsl:call-template name="dnl"/>
  </xsl:if>
</xsl:template>

<xsl:template match="usage">
  <xsl:text>&#10;&#10;</xsl:text>
  <xsl:text>\begin{tabular}{ll}&#10;</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>\end{tabular}&#10;&#10;</xsl:text>
</xsl:template>

<xsl:template match="arguments">
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
  <xsl:text> &amp; </xsl:text>
  <xsl:apply-templates/>
 <xsl:text>\\&#xa;</xsl:text>
</xsl:template>

<xsl:template match="argument">
  <xsl:if test="(@optional)">
    <xsl:text>\texttt{[ }</xsl:text>
  </xsl:if>
  <xsl:if test="(@separated)">
    <xsl:text>\texttt{; }</xsl:text>
  </xsl:if>
  <xsl:if test="(@alternate)">
    <xsl:call-template name="gettext">
      <xsl:with-param name="key" select="'or'"/>
    </xsl:call-template>
  </xsl:if>
  <xsl:if test="@flag">
    <xsl:text>\texttt{</xsl:text>
    <xsl:value-of select="@flag"/>
    <xsl:text>}</xsl:text>
  </xsl:if> 
  <xsl:text>\textsl{</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>} </xsl:text>
  <xsl:if test="(@optional)">
    <xsl:text>\texttt{ ]}</xsl:text>
  </xsl:if>
</xsl:template>

<xsl:template match="option">
  <xsl:if test="position() = 1">
    <xsl:choose>
      <xsl:when test="count(../option) > 1">
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
  </xsl:if>
  <xsl:text> &amp; </xsl:text>
  <xsl:apply-templates/>
  <xsl:text>&#10;</xsl:text>
</xsl:template>

<xsl:template match="example">
  <xsl:if test="position() = 1">
    <xsl:choose>
      <xsl:when test="count(../option) > 1">
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
  </xsl:if>
  <xsl:text> &amp; \texttt{</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>} \\ &#10;</xsl:text>
</xsl:template>

<xsl:template match="flag|argpunct">
  <xsl:text>\texttt{</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>}</xsl:text>
</xsl:template>

<xsl:template match="effect">
  <xsl:text> (</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>) \\</xsl:text>
</xsl:template>

<xsl:template match="demos">
  <xsl:if test="position() = 1">
    <xsl:choose>
      <xsl:when test="count(../example) > 1">
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
  </xsl:if>
  <xsl:if test="count(../example) > 1">
    <xsl:call-template name="gettext">
      <xsl:with-param name="key" select="'Seealso'"/>
    </xsl:call-template>
  </xsl:if>
  <xsl:text> &amp; </xsl:text>
  <xsl:apply-templates/>
  <xsl:text>&#10;</xsl:text>
</xsl:template>

<xsl:template match="demo">
  <xsl:if test="position() > 1">
    <xsl:text>, </xsl:text>
  </xsl:if>
  <xsl:text>\texttt{</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>}</xsl:text>
</xsl:template>

<xsl:template match="fnargs">
  <xsl:choose>
    <xsl:when test="count(fnarg) > 1">
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
  <xsl:text>\end{tabular}&#10;</xsl:text>
</xsl:template>

<xsl:template match="fnarg">
  <xsl:choose>
    <xsl:when test="position() > 1">
      <xsl:text>           &amp; </xsl:text>
    </xsl:when>
    <xsl:otherwise>
      <xsl:text>&amp; </xsl:text>
    </xsl:otherwise>
  </xsl:choose>
  <xsl:text>\textsl{</xsl:text>
  <xsl:if test="substring(., 1, 1) = '&amp;'">
    <xsl:text>\</xsl:text>
  </xsl:if>
  <xsl:apply-templates/>
  <xsl:text>}</xsl:text>
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
        <xsl:text>\texttt{null}</xsl:text>
      </xsl:when>
      <xsl:otherwise>
        <xsl:text>, </xsl:text>
        <xsl:call-template name="gettext">
          <xsl:with-param name="key" select="'optional'"/>
        </xsl:call-template>
      </xsl:otherwise>
    </xsl:choose>
  </xsl:if> 
  <xsl:text>)</xsl:text>
  <xsl:text>\\&#10;</xsl:text>
</xsl:template>

<xsl:template match="book">
  <xsl:text>\emph{</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>}</xsl:text>
</xsl:template>

<xsl:template match="quote">
  <xsl:text>``</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>''</xsl:text>
</xsl:template>

<xsl:template match="para">
  <xsl:if test="not(@context) or @context=$hlp or @context='cli'">
    <xsl:if test="@role='preformatted'">
      <xsl:text>\begin{raggedright}</xsl:text>
    </xsl:if>
    <xsl:apply-templates/>
    <xsl:if test="@role='preformatted'">
      <xsl:text>\end{raggedright}</xsl:text>
    </xsl:if>
    <xsl:if test="not(ancestor::entry)">
      <xsl:text>&#10;&#10;</xsl:text>
    </xsl:if>
  </xsl:if>
</xsl:template>

<xsl:template match="ilist">
  <xsl:if test="not(@context) or @context=$hlp">
    <xsl:text>&#10;\begin{itemize}&#10;</xsl:text>
    <xsl:apply-templates/>
    <xsl:text>\end{itemize}&#10;&#10;</xsl:text>
  </xsl:if>
</xsl:template>

<xsl:template match="olist">
  <xsl:text>&#10;\begin{enumerate}&#10;</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>\end{enumerate}&#10;&#10;</xsl:text>
</xsl:template>

<xsl:template match="li">
  <xsl:text>\item </xsl:text>
  <xsl:apply-templates/>
  <xsl:text>&#10;</xsl:text>
</xsl:template>

<xsl:template match="math">
  <xsl:text>$</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>$</xsl:text>
</xsl:template>

<xsl:template match="argname">
  <xsl:text>\textsl{</xsl:text>
  <xsl:value-of select="translate(., '_', '-')"/>
  <xsl:text>}</xsl:text>
</xsl:template>

<xsl:template match="super">
  <xsl:text>\ensuremath{^{</xsl:text>  
  <xsl:apply-templates/>
  <xsl:text>}}</xsl:text> 
</xsl:template>

<xsl:template match="sub">
  <xsl:text>\ensuremath{_{</xsl:text>  
  <xsl:apply-templates/>
  <xsl:text>}}</xsl:text> 
</xsl:template>

<xsl:template match="code">
  <xsl:text>&#10;\begin{code}</xsl:text>  
  <xsl:apply-templates/>
  <xsl:choose>
    <xsl:when test="substring(., string-length(.)-6, 1) = '&#10;'">
      <xsl:text>\end{code}&#10;</xsl:text>
    </xsl:when>
    <xsl:otherwise>
      <xsl:text>&#10;\end{code}&#10;</xsl:text> 
    </xsl:otherwise>
  </xsl:choose>
</xsl:template>

<xsl:template match="lit">
  <xsl:variable name="escape">
    <xsl:choose>
      <xsl:when test="contains(.,'_') or contains(.,'\') or 
                      contains(.,'$') or contains(.,'^') or
                      contains(.,'%') or contains(.,'&amp;') or
                      contains(.,'#') or contains(.,'{')">
        <xsl:text>yes</xsl:text>
      </xsl:when>    
      <xsl:otherwise>
        <xsl:text>no</xsl:text>
      </xsl:otherwise>
    </xsl:choose>    
  </xsl:variable>
  <xsl:choose>
    <xsl:when test="$escape='yes'">
      <xsl:text>\verb@</xsl:text>
    </xsl:when>    
    <xsl:otherwise>
      <xsl:text>\texttt{</xsl:text>
    </xsl:otherwise>
  </xsl:choose>
  <xsl:value-of select="translate(., '&#10;', '')"/>
  <xsl:choose>
    <xsl:when test="$escape='yes'">
      <xsl:text>@</xsl:text>
    </xsl:when>    
    <xsl:otherwise>
      <xsl:text>}</xsl:text>
    </xsl:otherwise>
  </xsl:choose>
</xsl:template>

<xsl:template match="equation">
  <xsl:value-of select="normalize-space(@tex)"/>
</xsl:template>

<xsl:template match="fncref">
  <xsl:text>\hyperlink{func-</xsl:text>
  <xsl:call-template name="escdollar">
    <xsl:with-param name="src" select="@targ"/>
  </xsl:call-template>
  <xsl:text>}</xsl:text>
  <xsl:text>{</xsl:text>
  <xsl:choose>
    <xsl:when test="@label">
      <xsl:value-of select="@label"/>
    </xsl:when>
    <xsl:otherwise>
      <xsl:call-template name="escdollar">
        <xsl:with-param name="src" select="@targ"/>
      </xsl:call-template>
    </xsl:otherwise>
  </xsl:choose>
  <xsl:text>}</xsl:text>
  <xsl:if test="parent::seelist and following-sibling::*">
    <xsl:text>, </xsl:text>
  </xsl:if>
</xsl:template>

<xsl:template match="guideref">
  <xsl:choose>
    <xsl:when test="$useGUG='true'">
      <xsl:text>\GUG{}</xsl:text>
    </xsl:when>
    <xsl:otherwise>
      <xref linkend="{@targ}"/>
    </xsl:otherwise>
  </xsl:choose>
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

<xsl:template match="menu-path">
  <xsl:call-template name="gettext">
    <xsl:with-param name="key" select="'menupath'"/>
  </xsl:call-template>
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="description">
  <xsl:apply-templates/>
</xsl:template>

</xsl:stylesheet>
