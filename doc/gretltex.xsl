<?xml version="1.0"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" >

<xsl:param name="lang" select="'en'"/>

<xsl:output method="text" encoding="iso-8859-1"/>

<xsl:include href="xsl/common.xsl"/>
<xsl:include href="xsl/common/common.xsl"/>
<xsl:include href="xsl/common/l10n.xsl"/>
<xsl:include href="xsl/xref.xsl"/>

<xsl:strip-space elements="title"/>
<xsl:preserve-space elements="*"/>

<xsl:template match="chapter">
  <xsl:text>\chapter{</xsl:text>
  <xsl:value-of select="title"/>
  <xsl:text>}&#10;\label{</xsl:text>
  <xsl:value-of select="@id"/>
  <xsl:text>}&#10;&#10;</xsl:text>
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="appendix">
  <xsl:text>\chapter{</xsl:text>
  <xsl:value-of select="title"/>
  <xsl:text>}&#10;\label{</xsl:text>
  <xsl:value-of select="@id"/>
  <xsl:text>}&#10;&#10;</xsl:text>
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="sect1">
  <xsl:text>&#10;&#10;\section{</xsl:text>
  <xsl:value-of select="title"/>
  <xsl:text>}&#10;</xsl:text>
  <xsl:if test="@id">
    <xsl:text>\label{</xsl:text>
    <xsl:value-of select="@id"/>
    <xsl:text>}</xsl:text>
  </xsl:if>
  <xsl:text>&#10;&#10;</xsl:text>
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="title"/>

<xsl:template match="sect2">
  <xsl:text>&#10;&#10;\subsection{</xsl:text>
  <xsl:value-of select="title"/>
  <xsl:text>}&#10;</xsl:text>
  <xsl:choose>
    <xsl:when test="@xreflabel">
      <xsl:text>\hypertarget{</xsl:text>
      <xsl:value-of select="@id"/>
      <xsl:text>}{}</xsl:text>    
    </xsl:when>
    <xsl:otherwise>
      <xsl:if test="@id">
        <xsl:text>\label{</xsl:text>
        <xsl:value-of select="@id"/>
        <xsl:text>}</xsl:text>
      </xsl:if>    
    </xsl:otherwise>    
  </xsl:choose>
  <xsl:text>&#10;&#10;</xsl:text>
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="command">
  <xsl:text>\cmd{</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>}</xsl:text>
</xsl:template>

<xsl:template match="emphasis">
  <xsl:text>\emph{</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>}</xsl:text>
</xsl:template>

<xsl:template match="function">
  <xsl:text>\cmd{</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>}</xsl:text>
</xsl:template>

<xsl:template match="citetitle">
  <xsl:text>\emph{</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>}</xsl:text>
</xsl:template>

<xsl:template match="application">
  <xsl:text>\app{</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>}</xsl:text>
</xsl:template>

<xsl:template match="quote">
  <xsl:text>``</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>''</xsl:text>
</xsl:template>

<xsl:template match="para">
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
</xsl:template>

<xsl:template match="variablelist">
  <xsl:text>&#10;\begin{description}&#10;</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>\end{description}&#10;&#10;</xsl:text>
</xsl:template>

<xsl:template match="varlistentry">
  <xsl:text>\item[</xsl:text>
  <xsl:value-of select="term"/>
  <xsl:text>] </xsl:text>
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="term"/>

<xsl:template match="itemizedlist">
  <xsl:text>&#10;\begin{itemize}&#10;</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>\end{itemize}&#10;&#10;</xsl:text>
</xsl:template>

<xsl:template match="orderedlist">
  <xsl:text>&#10;\begin{enumerate}&#10;</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>\end{enumerate}&#10;&#10;</xsl:text>
</xsl:template>

<xsl:template match="listitem">
  <xsl:if test="not(parent::varlistentry)">
    <xsl:text>\item </xsl:text>
  </xsl:if>
  <xsl:apply-templates/>
  <xsl:text>&#10;</xsl:text>
</xsl:template>

<xsl:template match="mathvar">
  <xsl:text>$</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>$</xsl:text>
</xsl:template>

<xsl:template match="book">
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="replaceable">
  <xsl:text>\textsl{</xsl:text>
  <xsl:value-of select="translate(., '_', '-')"/>
  <xsl:text>}</xsl:text>
</xsl:template>

<xsl:template match="footnote">
  <xsl:text>\footnote{</xsl:text>  
  <xsl:apply-templates/>
  <xsl:text>}</xsl:text> 
</xsl:template>

<xsl:template match="superscript">
  <xsl:text>\ensuremath{^{</xsl:text>  
  <xsl:apply-templates/>
  <xsl:text>}}</xsl:text> 
</xsl:template>

<xsl:template match="subscript">
  <xsl:text>\ensuremath{_{</xsl:text>  
  <xsl:apply-templates/>
  <xsl:text>}}</xsl:text> 
</xsl:template>

<xsl:template match="programlisting">
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

<xsl:template match="literal|filename|varname">
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

<xsl:template match="guimenu|guimenuitem">
  <xsl:text>\textsf{</xsl:text>  
  <xsl:apply-templates/>
  <xsl:text>}</xsl:text> 
</xsl:template>

<xsl:template match="informalequation">
  <xsl:value-of select="alt[@role='tex']"/>
</xsl:template>

<xsl:template match="inlineequation">
  <xsl:value-of select="alt[@role='tex']"/>
</xsl:template>

<xsl:template match="figure">
  <xsl:text>&#10;\begin{figure}[htbp]&#10;</xsl:text>
  <xsl:text>\begin{center}&#10;</xsl:text>
  <xsl:text>\includegraphics[scale=0.5]{</xsl:text>
  <xsl:choose>
    <xsl:when test="screenshot">
      <xsl:apply-templates select="screenshot"/>
    </xsl:when> 
    <xsl:when test="graphic">
      <xsl:apply-templates select="graphic"/>
    </xsl:when> 
  </xsl:choose>
  <xsl:text>}&#10;</xsl:text>
  <xsl:text>\end{center}&#10;</xsl:text>
  <xsl:text>\caption{</xsl:text>
  <xsl:value-of select="title"/>
  <xsl:text>}&#10;</xsl:text>  
  <xsl:text>\label{</xsl:text>
  <xsl:value-of select="@id"/>
  <xsl:text>}&#10;</xsl:text>
  <xsl:text>\end{figure}&#10;</xsl:text>
</xsl:template>

<xsl:template match="screenshot">
  <xsl:choose>
    <xsl:when test="graphic">
      <xsl:apply-templates select="graphic"/>
    </xsl:when>
  </xsl:choose> 
</xsl:template>

<xsl:template match="graphic">
  <xsl:if test="not(ancestor::informalequation) and not(ancestor::equation)">
    <xsl:choose>
      <xsl:when test="ancestor::screenshot or ancestor::figure">
        <xsl:value-of select="@fileref"/>
      </xsl:when>
      <xsl:otherwise>
        <xsl:text>\begin{center}&#10;</xsl:text>
        <xsl:text>\includegraphics[scale=0.5]{</xsl:text>
        <xsl:value-of select="@fileref"/>
        <xsl:text>}&#10;\end{center}&#10;</xsl:text>
      </xsl:otherwise>
    </xsl:choose> 
  </xsl:if>
</xsl:template>

<xsl:template match="tip">
  <xsl:text>&#10;\tip{</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>}&#10;</xsl:text>
</xsl:template>

<xsl:template match="example">
  <xsl:text>&#10;&#10;\begin{script}[htbp]&#10;</xsl:text>
  <xsl:text>\caption{</xsl:text>
  <xsl:value-of select="title"/>
  <xsl:text>}&#10;</xsl:text>
  <xsl:text>\label{</xsl:text>
  <xsl:value-of select="@id"/>
  <xsl:text>}&#10;</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>\end{script}&#10;&#10;</xsl:text>
</xsl:template>

<xsl:template match="entry">
  <xsl:choose>
    <xsl:when test="ancestor::thead">
      <xsl:text>\textit{</xsl:text>
      <xsl:apply-templates/>
      <xsl:text>}</xsl:text>
    </xsl:when>
    <xsl:otherwise>
      <xsl:apply-templates/>
    </xsl:otherwise>
  </xsl:choose>
  <xsl:choose>
    <xsl:when test="following-sibling::entry">
      <xsl:text> &amp; </xsl:text>
    </xsl:when>
  </xsl:choose>
</xsl:template>

<xsl:template match="row">
  <xsl:apply-templates/>
  <xsl:choose>
    <xsl:when test="ancestor::thead">
      <xsl:text> \\ [4pt]&#10;</xsl:text>
    </xsl:when>
    <xsl:otherwise>
      <xsl:text> \\&#10;</xsl:text>
    </xsl:otherwise>
  </xsl:choose>
</xsl:template>

<xsl:template match="tbody">
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="thead">
  <xsl:apply-templates/>
</xsl:template>

<xsl:template name="tabcolstring">
  <xsl:param name="colcount"/>  
  <xsl:call-template name="tab-col-loop">
    <xsl:with-param name="i" select="1"/>
    <xsl:with-param name="toid" select="$colcount"/>
    <xsl:with-param name="stepid" select="1"/>
    <xsl:with-param name="colcount" select="$colcount"/>
  </xsl:call-template>
</xsl:template>

<xsl:template match="colspec">
  <xsl:choose>
    <xsl:when test="parent::*[@style='allpara']">
      <xsl:text>&gt;{\raggedright\arraybackslash}p{</xsl:text>
      <xsl:value-of select="@colwidth"/>
      <xsl:text>}&#10;</xsl:text>   
      <xsl:if test="not(following-sibling::colspec)">
        <xsl:text>}&#10;</xsl:text>
      </xsl:if>  
    </xsl:when>
    <xsl:when test="parent::*[@style='lpara']">
      <xsl:choose>
        <xsl:when test="position() = 1">
          <xsl:text>&gt;{\raggedright\arraybackslash}p{</xsl:text>
          <xsl:value-of select="@colwidth"/>
          <xsl:text>}&#10;</xsl:text>
        </xsl:when>
        <xsl:otherwise>
          <xsl:text>l</xsl:text>
        </xsl:otherwise>
      </xsl:choose>
      <xsl:if test="not(following-sibling::colspec)">
        <xsl:text>}&#10;</xsl:text>
      </xsl:if>
    </xsl:when>
    <xsl:when test="parent::*[@style='rpara']">
      <xsl:choose>
        <xsl:when test="position() = 2">
          <xsl:text>&gt;{\raggedright\arraybackslash}p{</xsl:text>
          <xsl:value-of select="@colwidth"/>
          <xsl:text>}&#10;</xsl:text>
        </xsl:when>
        <xsl:otherwise>
          <xsl:text>l</xsl:text>
        </xsl:otherwise>
      </xsl:choose>
      <xsl:if test="not(following-sibling::colspec)">
        <xsl:text>}&#10;</xsl:text>
      </xsl:if>
    </xsl:when>
  </xsl:choose>
</xsl:template>

<xsl:template match="tgroup">
  <xsl:text>\begin{tabular}{</xsl:text>
  <xsl:if test="not(@style='lpara' or @style='rpara' or @style='allpara')">
    <xsl:call-template name="tabcolstring">
      <xsl:with-param name="colcount" select="@cols"/>
    </xsl:call-template>      
    <xsl:text>}&#10;</xsl:text>
  </xsl:if>
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="informaltable">
  <xsl:text>&#10;&#10;</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>\end{tabular}&#10;&#10;</xsl:text>
</xsl:template>

<xsl:template match="table">
  <xsl:text>&#10;&#10;\begin{table}[htbp]&#10;</xsl:text>
  <xsl:text>\caption{</xsl:text>
  <xsl:value-of select="title"/>
  <xsl:text>}&#10;</xsl:text>
  <xsl:text>\label{</xsl:text>
  <xsl:value-of select="@id"/>
  <xsl:text>}&#10;</xsl:text>  
  <xsl:text>\begin{center}&#10;</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>\end{tabular}&#10;</xsl:text>
  <xsl:text>\end{center}&#10;</xsl:text>
  <xsl:text>\end{table}&#10;&#10;</xsl:text>
</xsl:template>

<xsl:template name="tab-col-loop">
  <xsl:param name="i"/>
  <xsl:param name="toid"/>
  <xsl:param name="stepid"/>
  <xsl:param name="colcount"/>
    <xsl:text>l</xsl:text>
  <xsl:if test="$i+$stepid &lt;= $toid">
    <xsl:call-template name="tab-col-loop">
      <xsl:with-param name="i" select="$i + $stepid"/>
      <xsl:with-param name="toid" select="$toid"/>
      <xsl:with-param name="stepid" select="$stepid"/>
      <xsl:with-param name="colcount" select="$colcount"/>
    </xsl:call-template>
  </xsl:if>
</xsl:template>

<xsl:template match="bibliomset/title">
  <xsl:variable name="relation" select="../@relation"/>
  <xsl:choose>
    <xsl:when test="$relation='article'">
      <xsl:text>``</xsl:text>
      <xsl:apply-templates/>
      <xsl:text>''</xsl:text>
    </xsl:when>
    <xsl:otherwise>
      <xsl:text>\emph{</xsl:text>
      <xsl:apply-templates/>
      <xsl:text>}</xsl:text>
    </xsl:otherwise>
  </xsl:choose>
</xsl:template>

<xsl:template match="bibliomset">
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="bibliomixed/title">
  <xsl:text>\emph{</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>}</xsl:text>
</xsl:template>

<xsl:template match="bibliomixed">
  <xsl:apply-templates/>
  <xsl:text>&#10;</xsl:text>
</xsl:template>

<xsl:template match="bibliography">
  <xsl:text>\begin{thebibliography}&#10;</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>&#10;\end{thebibliography}&#10;</xsl:text>
</xsl:template>

<xsl:template match="bookinfo"/>

</xsl:stylesheet>
