<!-- DSSSL Stylesheet fragment TeXMath.dsl -->

(element (informalequation texmath)
 (make display-group
   space-before: (* %para-sep% 0.0)
   space-after: %block-sep%
   start-indent: %body-start-indent%
   (literal "BEGINTEXMATH")
   (literal (data (current-node)))
   (literal "ENDTEXMATH")))

(element (inlineequation texmath)
 (make sequence
   (literal "BEGINTEXMATH")
   (literal (data (current-node)))
   (literal "ENDTEXMATH")))   

(element (informalequation graphic)
  (empty-sosofo))

(element (inlineequation inlinegraphic)
  (empty-sosofo))
  

