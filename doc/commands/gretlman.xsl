<?xml version="1.0"?>
<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform" >

<xsl:template match="commandlist"> 
 <xsl:apply-templates/> 
</xsl:template>

<xsl:template match="command">
  <xsl:text>&#xa;</xsl:text>
  <sect2 id="cmd-{@name}" xreflabel="{@name}">
    <title><xsl:value-of select="@name"/></title>
  <xsl:text>&#xa;</xsl:text>
  <xsl:apply-templates/>
  <xsl:text>&#xa;</xsl:text>
  </sect2>
  <xsl:text>&#xa;</xsl:text>
</xsl:template>

<xsl:template match="usage">
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

<xsl:template match="arguments">
 <row>
   <xsl:choose>
     <xsl:when test="count(argument) > 1">
       <entry>Arguments:</entry>
     </xsl:when>
     <xsl:otherwise>
       <entry>Argument:</entry>
     </xsl:otherwise>
   </xsl:choose>
  <entry><xsl:apply-templates/></entry>
 </row> 
 <xsl:text>&#xa;</xsl:text>
</xsl:template>

<xsl:template match="argument">
<xsl:if test="(@separated)">; </xsl:if>
<xsl:if test="(@alternate)"> or </xsl:if>
<xsl:if test="(@optional)"><literal>[ </literal></xsl:if>
  <xsl:if test="@flag">
    <literal><xsl:value-of select="@flag"/></literal>
  </xsl:if> 
  <replaceable><xsl:apply-templates/></replaceable>
  <xsl:text> </xsl:text>
<xsl:if test="(@optional)"><literal>] </literal></xsl:if>
</xsl:template>

<xsl:template match="option">
  <row>
    <xsl:choose>
      <xsl:when test="position() = 1">
        <xsl:choose>
          <xsl:when test="count(../option) > 1">
            <entry>Options:</entry>
          </xsl:when>
          <xsl:otherwise>
            <entry>Option:</entry>
          </xsl:otherwise>
        </xsl:choose>
      </xsl:when>
      <xsl:otherwise>
        <entry/>
      </xsl:otherwise>
    </xsl:choose>
    <entry><xsl:apply-templates/></entry>
  </row>
  <xsl:text>&#xa;</xsl:text>
</xsl:template>

<xsl:template match="example">
  <row>
    <xsl:choose>
      <xsl:when test="position() = 1">
        <xsl:choose>
          <xsl:when test="count(../example) > 1">
            <entry>Examples:</entry>
          </xsl:when>
          <xsl:otherwise>
            <entry>Example:</entry>
          </xsl:otherwise>
        </xsl:choose>
      </xsl:when>
      <xsl:otherwise>
        <entry/>
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

<xsl:template match="super">
  <superscript><xsl:apply-templates/></superscript>
</xsl:template>

<xsl:template match="sub">
  <subscript><xsl:apply-templates/></subscript>
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

<xsl:template match="ilist">
  <itemizedlist><xsl:apply-templates/></itemizedlist>
  <xsl:text>&#xa;</xsl:text>
</xsl:template>

<xsl:template match="li">
  <listitem><xsl:apply-templates/></listitem>
  <xsl:text>&#xa;</xsl:text>
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

<xsl:template match="cmdref">
  <xref linkend="cmd-{@targ}"/>
</xsl:template>

<xsl:template match="manref">
<xsl:value-of select="@before"/>
<xsl:text> </xsl:text><xref linkend="{@targ}"/>
<xsl:text> </xsl:text><xsl:value-of select="@after"/>
</xsl:template>

<xsl:template match="tabref">
  <xref linkend="{@targ}"/>
</xsl:template>

<xsl:template match="menu-path">
  <para>Menu path: <xsl:apply-templates/></para>
</xsl:template>

<xsl:template match="other-access">
  <para>Other access: <xsl:apply-templates/></para>
</xsl:template>

</xsl:stylesheet>
