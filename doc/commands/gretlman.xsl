<?xml version="1.0"?>
<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform" >

<xsl:param name="hlp">cli</xsl:param>
<xsl:param name="standalone">false</xsl:param>
<xsl:param name="lang" select="'en'"/>

<xsl:output method="xml" omit-xml-declaration="yes" indent="yes"/>

<xsl:variable name="phrases"
  select="document('manual_strings.xml')/phrases"/>

<xsl:template name="gettext">
  <xsl:param name="key"/>
  <xsl:value-of
    select="$phrases/phrase[@key=$key and @lang=$lang]"/>
</xsl:template>

<xsl:template match="commandlist"> 
 <xsl:apply-templates/> 
</xsl:template>

<xsl:template match="command[not(@context) or @context=$hlp]">
  <xsl:text>&#xa;</xsl:text>
  <sect2 id="cmd-{@name}" xreflabel="{@name}">
    <title><xsl:value-of select="@name"/></title>
    <xsl:call-template name="nl"/>
    <xsl:apply-templates/>
    <xsl:call-template name="nl"/>
  </sect2>
  <xsl:call-template name="nl"/>
</xsl:template>

<xsl:template match="command[@context and @context!=$hlp]"/>

<xsl:template match="description[not(@context) or @context=$hlp]">
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="description[@context and @context!=$hlp]"/>

<xsl:template match="usage">
  <xsl:text>&#xa;</xsl:text>
  <informaltable role="cmd" frame="none">
    <tgroup cols="2"><colspec colnum="1" colwidth="82pt"/>
    <tbody>
      <xsl:call-template name="nl"/>
      <xsl:apply-templates/>
    </tbody>
  </tgroup>
</informaltable> 
<xsl:call-template name="nl"/>
</xsl:template>

<xsl:template match="arguments">
 <row>
   <xsl:choose>
     <xsl:when test="count(argument) > 1">
       <entry>
         <xsl:call-template name="gettext">
           <xsl:with-param name="key" select="'args'"/>
         </xsl:call-template>
       </entry>
     </xsl:when>
     <xsl:otherwise>
       <entry>
         <xsl:call-template name="gettext">
           <xsl:with-param name="key" select="'arg'"/>
         </xsl:call-template>
       </entry>
     </xsl:otherwise>
   </xsl:choose>
  <entry><xsl:apply-templates/></entry>
 </row> 
 <xsl:text>&#xa;</xsl:text>
</xsl:template>

<xsl:template match="argument">
<xsl:if test="(@separated)">; </xsl:if>
<xsl:if test="(@alternate)">
  <xsl:call-template name="gettext">
    <xsl:with-param name="key" select="'or'"/>
  </xsl:call-template>
</xsl:if>
<xsl:if test="(@optional)"><literal>[ </literal></xsl:if>
  <xsl:if test="@flag">
    <literal><xsl:value-of select="@flag"/></literal>
  </xsl:if> 
  <replaceable><xsl:apply-templates/></replaceable>
  <xsl:text> </xsl:text>
<xsl:if test="(@optional)"><literal>] </literal></xsl:if>
</xsl:template>

<xsl:template match="altform">
  <row>
    <xsl:choose>
      <xsl:when test="position() = 1">
        <entry>
          <xsl:call-template name="gettext">
            <xsl:with-param name="key" select="'altforms'"/>
          </xsl:call-template>
        </entry>
      </xsl:when>
      <xsl:otherwise>
        <entry/>
      </xsl:otherwise>
    </xsl:choose>
    <entry><xsl:apply-templates/></entry>
  </row>
  <xsl:text>&#xa;</xsl:text>
</xsl:template>

<xsl:template match="option">
  <row>
    <xsl:choose>
      <xsl:when test="position() = 1">
        <xsl:choose>
          <xsl:when test="count(../option) > 1">
            <entry>
              <xsl:call-template name="gettext">
                <xsl:with-param name="key" select="'opts'"/>
              </xsl:call-template>
            </entry>
          </xsl:when>
          <xsl:otherwise>
            <entry>
              <xsl:call-template name="gettext">
                <xsl:with-param name="key" select="'opt'"/>
              </xsl:call-template>
            </entry>
          </xsl:otherwise>
        </xsl:choose>
      </xsl:when>
      <xsl:otherwise>
        <entry/>
      </xsl:otherwise>
    </xsl:choose>
    <entry><xsl:apply-templates/></entry>
  </row>
  <xsl:call-template name="nl"/>
</xsl:template>

<xsl:template match="example">
  <row>
    <xsl:choose>
      <xsl:when test="position() = 1">
        <xsl:choose>
          <xsl:when test="count(../example) > 1">
            <entry>
              <xsl:call-template name="gettext">
                <xsl:with-param name="key" select="'examples'"/>
              </xsl:call-template>
            </entry>
          </xsl:when>
          <xsl:otherwise>
            <entry>
              <xsl:call-template name="gettext">
                <xsl:with-param name="key" select="'example'"/>
              </xsl:call-template>
            </entry>
          </xsl:otherwise>
        </xsl:choose>
      </xsl:when>
      <xsl:otherwise>
        <entry/>
      </xsl:otherwise>
    </xsl:choose>
    <entry><literal><xsl:apply-templates/></literal></entry>
  </row>
  <xsl:call-template name="nl"/>
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

<xsl:template match="lit">
  <literal><xsl:apply-templates/></literal>
</xsl:template>

<xsl:template match="mathvar">
  <emphasis><xsl:apply-templates/></emphasis>
</xsl:template>

<xsl:template match="emphasis">
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

<xsl:template match="footnote">
  <footnote><xsl:apply-templates/></footnote>
</xsl:template>

<xsl:template match="filename">
  <filename><xsl:apply-templates/></filename>
</xsl:template>

<xsl:template match="func">
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

<xsl:template match="ilist[not(@context) or @context=$hlp]">
  <itemizedlist><xsl:apply-templates/></itemizedlist>
  <xsl:call-template name="nl"/>
</xsl:template>

<xsl:template match="ilist[@context and @context!=$hlp]"/>

<xsl:template match="li">
  <listitem><xsl:apply-templates/></listitem>
  <xsl:call-template name="nl"/>
</xsl:template>

<xsl:template match="para[not(@context) or @context=$hlp]">
  <xsl:call-template name="nl"/>
  <para>
    <xsl:apply-templates/>
  </para>
  <xsl:call-template name="nl"/>
</xsl:template>

<xsl:template match="para[@context and @context!=$hlp]"/>

<xsl:template match="code">
<programlisting>
<xsl:apply-templates/></programlisting>
</xsl:template>

<xsl:template match="cmdref">
  <xref linkend="cmd-{@targ}"/>
</xsl:template>

<xsl:template match="manref">
  <xsl:choose>
    <xsl:when test="$standalone='true'">
      <xsl:call-template name="gettext">
        <xsl:with-param name="key" select="'fullman'"/>
      </xsl:call-template>
    </xsl:when>
    <xsl:otherwise>
      <xref linkend="{@targ}"/>
    </xsl:otherwise>
  </xsl:choose>
</xsl:template>

<xsl:template match="tabref">
  <xref linkend="{@targ}"/>
</xsl:template>

<xsl:template match="menu-path">
  <para>
    <xsl:call-template name="gettext">
      <xsl:with-param name="key" select="'menupath'"/>
    </xsl:call-template>
  <xsl:apply-templates/></para>
</xsl:template>

<xsl:template match="other-access">
  <para>
    <xsl:call-template name="gettext">
      <xsl:with-param name="key" select="'otheraccess'"/>
    </xsl:call-template>
  <xsl:apply-templates/></para>
</xsl:template>

<xsl:template match="table[not(@context) or @context=$hlp]">
  <xsl:choose>
    <xsl:when test="@id">
      <table id="{@id}" frame="none">
        <title><xsl:value-of select="@title"/></title>
        <tgroup cols="2">
          <thead>
            <row>
              <entry><xsl:value-of select="@lhead"/></entry>
              <entry><xsl:value-of select="@rhead"/></entry>
            </row>
          </thead>
          <tbody>
            <xsl:apply-templates/>
          </tbody>
        </tgroup>
      </table>
    </xsl:when>
    <xsl:otherwise>
      <informaltable role="cmd" frame="none">
        <tgroup cols="2">
          <xsl:if test="@lwidth and @rwidth">
            <colspec colwidth="{@lwidth}"/>
            <colspec colwidth="{@rwidth}"/>
          </xsl:if>
          <tbody>
            <xsl:apply-templates/>
          </tbody>
        </tgroup>
      </informaltable>
    </xsl:otherwise>
  </xsl:choose>
</xsl:template>

<xsl:template match="table[@context and @context!=$hlp]"/>

<xsl:template match="row">
  <row><xsl:apply-templates/></row>
</xsl:template>

<xsl:template match="cell">
  <entry><xsl:apply-templates/></entry>
</xsl:template>

<xsl:template name="nl">
  <xsl:text>&#10;</xsl:text>  
</xsl:template>

</xsl:stylesheet>
