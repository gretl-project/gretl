<?xml version="1.0"?>
<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform" >

<xsl:template match="commandlist"> 
 <xsl:apply-templates/> 
</xsl:template>

<xsl:template match="command">
  <xsl:text>&#xa;</xsl:text>
  <sect2 id="cmd-{@name}">
    <title><xsl:value-of select="@name"/></title>
  <xsl:text>&#xa;</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>&#xa;</xsl:text>
  </sect2>
  <xsl:text>&#xa;</xsl:text>
</xsl:template>

<xsl:template match="toptable">
 <xsl:text>&#xa;</xsl:text>
 <informaltable role="cmd" frame="none">
 <tgroup cols="2"><colspec colnum="1" colwidth="82pt"/>
 <tbody>
  <xsl:text>&#xa;</xsl:text>
    <xsl:apply-templates/>
    </tbody>
   </tgroup>
  </informaltable> 
 <xsl:text>&#xa;</xsl:text>
</xsl:template>

<xsl:template match="option">
 <row>
   <xsl:choose>
     <xsl:when test="position() = 1">
        <entry>Options:</entry>
     </xsl:when>
     <xsl:otherwise>
        <entry/>
     </xsl:otherwise>
    </xsl:choose>
    <entry><xsl:apply-templates/></entry>
  </row>
  <xsl:text>&#xa;</xsl:text>
</xsl:template>

<xsl:template match="arguments">
  <row>
    <entry>Arguments:</entry>
    <entry><replaceable><xsl:apply-templates/></replaceable></entry>
  </row>
  <xsl:text>&#xa;</xsl:text>
</xsl:template>

<xsl:template match="argument">
  <xsl:apply-templates/><xsl:text> </xsl:text>
</xsl:template>

<xsl:template match="example">
 <row>
 <xsl:choose>
   <xsl:when test="position() = 1">
     <entry>Examples:</entry>
    </xsl:when>
    <xsl:otherwise>
      <entry></entry>
    </xsl:otherwise>
 </xsl:choose>
 <entry><literal><xsl:apply-templates/></literal></entry>
 </row>
 <xsl:text>&#xa;</xsl:text>
</xsl:template>

<xsl:template match="flag">
  <literal><xsl:apply-templates/></literal>
</xsl:template>

<xsl:template match="effect">
  <xsl:text> (</xsl:text><xsl:apply-templates/><xsl:text>)</xsl:text>
</xsl:template>

<xsl:template match="repl">
  <replaceable><xsl:apply-templates/></replaceable>
</xsl:template>

<xsl:template match="cmd">
  <command><xsl:apply-templates/></command>
</xsl:template>

<xsl:template match="program">
  <application><xsl:apply-templates/></application>
</xsl:template>

<xsl:template match="literal">
  <literal><xsl:apply-templates/></literal>
</xsl:template>

<xsl:template match="mathvar">
  <emphasis><xsl:apply-templates/></emphasis>
</xsl:template>

<xsl:template match="book">
  <citetitle><xsl:apply-templates/></citetitle>
</xsl:template>

<xsl:template match="quote">
  <quote><xsl:apply-templates/></quote>
</xsl:template>

<xsl:template match="filename">
  <filename><xsl:apply-templates/></filename>
</xsl:template>

<xsl:template match="function">
  <function><xsl:apply-templates/></function>
</xsl:template>

<xsl:template match="equation">
 <xsl:choose>
  <xsl:when test="(@status='display')">
   <informalequation>
    <alt role="tex"><xsl:value-of select="@tex"/></alt>
    <graphic fileref="figures/{@graphic}"/>
   </informalequation>
  </xsl:when>
  <xsl:otherwise>
   <inlineequation>
    <alt role="tex"><xsl:value-of select="@tex"/></alt> 
    <inlinemediaobject>
     <imageobject>
      <imagedata align="center" fileref="figures/{@graphic}"/>
     </imageobject>
    </inlinemediaobject>
   </inlineequation>
  </xsl:otherwise>
 </xsl:choose>
</xsl:template>

<xsl:template match="para">
  <xsl:text>&#xa;</xsl:text>
  <para>
    <xsl:apply-templates/>
  </para>
  <xsl:text>&#xa;</xsl:text>
</xsl:template>

<xsl:template match="code">
<programlisting>
<xsl:apply-templates/></programlisting>
</xsl:template>

<xsl:template match="menu-path">
 <para>Menu path: <xsl:apply-templates/></para>
</xsl:template>

<xsl:template match="other-access">
  <para>Other access: <xsl:apply-templates/></para>
</xsl:template>

</xsl:stylesheet>
