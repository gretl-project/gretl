<!-- DSSSL Stylesheet fragment TeXMath.dsl 
     Allin Cottrell <cottrell@wfu.edu>, November 2001 -->

(element (equation graphic) (empty-sosofo))
(element (equation alt)
 (make display-group
   (literal "BEGINTEXLITERAL")
   (literal (data (current-node)))
   (literal "ENDTEXLITERAL")))

(element (informalequation graphic) (empty-sosofo))
(element (informalequation alt)
 (make display-group
   (literal "BEGINTEXLITERAL")
   (literal (data (current-node)))
   (literal "ENDTEXLITERAL")))
   
(element (inlineequation graphic) (empty-sosofo))
(element (inlineequation alt)
 (make sequence
   (literal "BEGINTEXLITERAL")
   (literal (data (current-node)))
   (literal "ENDTEXLITERAL")))


