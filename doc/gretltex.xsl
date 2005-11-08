<?xml version="1.0"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" >

<xsl:param name="lang" select="'en'"/>

<xsl:output method="text" encoding="iso-8859-1"/>

<xsl:template match="chapter">
  <xsl:text>\chapter{</xsl:text>
  <xsl:value-of select="title"/>
  <xsl:text>}&#10;\label{</xsl:text>
  <xsl:value-of select="@id"/>
  <xsl:text>}</xsl:text>
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="sect1">
  <xsl:text>\section{</xsl:text>
  <xsl:value-of select="title"/>
  <xsl:text>}</xsl:text>
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="title"/>

<xsl:template match="sect2">
  <xsl:text>\subsection{</xsl:text>
  <xsl:value-of select="title"/>
  <xsl:text>}</xsl:text>
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

<xsl:template match="variablelist">
  <xsl:text>\begin{description}&#10;</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>\end{description}&#10;</xsl:text>
</xsl:template>

<xsl:template match="varlistentry">
  <xsl:text>\item[</xsl:text>
  <xsl:value-of select="term"/>
  <xsl:text>] </xsl:text>
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="term"/>

<xsl:template match="itemizedlist">
  <xsl:text>\begin{itemize}&#10;</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>\end{itemize}&#10;</xsl:text>
</xsl:template>

<xsl:template match="orderedlist">
  <xsl:text>\begin{enumerate}&#10;</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>\end{enumerate}&#10;</xsl:text>
</xsl:template>

<xsl:template match="listitem">
  <xsl:if test="not(parent::varlistentry)">
    <xsl:text>\item </xsl:text>
  </xsl:if>
  <xsl:apply-templates/>
  <xsl:text>&#10;</xsl:text>
</xsl:template>

<xsl:template match="xref">
  <xsl:text>\ref{</xsl:text>
  <xsl:value-of select="@linkend"/>
  <xsl:text>}</xsl:text>
</xsl:template>

<xsl:template match="mathvar">
  <xsl:text>$</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>$</xsl:text>
</xsl:template>

<xsl:template match="book">
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="footnote">
  <xsl:text>\footnote{</xsl:text>  
  <xsl:apply-templates/>
  <xsl:text>}</xsl:text> 
</xsl:template>

<xsl:template match="superscript">
  <xsl:text>^{</xsl:text>  
  <xsl:apply-templates/>
  <xsl:text>}</xsl:text> 
</xsl:template>

<xsl:template match="subscript">
  <xsl:text>_{</xsl:text>  
  <xsl:apply-templates/>
  <xsl:text>}</xsl:text> 
</xsl:template>

<xsl:template match="programlisting">
  <xsl:text>\begin{verbatim}</xsl:text>  
  <xsl:apply-templates/>
  <xsl:text>\end{verbatim}</xsl:text> 
</xsl:template>

<xsl:template match="literal">
  <xsl:text>\verb+</xsl:text>  
  <xsl:apply-templates/>
  <xsl:text>+</xsl:text> 
</xsl:template>

<xsl:template match="filename">
  <xsl:text>\url{</xsl:text>  
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
  <xsl:text>\begin{figure}[htbp]&#10;</xsl:text>
  <xsl:text>\caption{</xsl:text>
  <xsl:value-of select="title"/>
  <xsl:text>}&#10;</xsl:text>  
  <xsl:text>\label{</xsl:text>
  <xsl:value-of select="@id"/>
  <xsl:text>}&#10;</xsl:text>
  <xsl:text>\centering&#10;</xsl:text>
  <xsl:text>\includegraphics[scale=0.5]{</xsl:text>
  <xsl:choose>
    <xsl:when test="screenshot">
      <xsl:apply-templates select="screenshot"/>
    </xsl:when>    
  </xsl:choose>
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
  <xsl:value-of select="@fileref"/>
</xsl:template>

</xsl:stylesheet>
