<!ELEMENT mathvar (#PCDATA|subscript|superscript)* >
<!-- Allow <mathvar> in paragraphs -->
<!ENTITY % local.para.char.mix
        "|mathvar">
