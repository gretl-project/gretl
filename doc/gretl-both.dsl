<!DOCTYPE style-sheet PUBLIC "-//James Clark//DTD DSSSL Style Sheet//EN" [
<!ENTITY % html "IGNORE">
<![%html;[
<!ENTITY % print "IGNORE">
<!ENTITY docbook.dsl SYSTEM "/usr/share/sgml/docbook-dsssl-1.73/html/docbook.dsl" CDATA dsssl>
]]>
<!ENTITY % print "INCLUDE">
<![%print;[
<!ENTITY docbook.dsl SYSTEM "/usr/share/sgml/docbook-dsssl-1.73/print/docbook.dsl" CDATA dsssl>
]]>
<!ENTITY htmlmath.dsl SYSTEM "HTMLMath.dsl">
<!ENTITY texmath.dsl SYSTEM "TeXMath.dsl">
<!ENTITY titlepage.dsl SYSTEM "titlepage.dsl">
]>

<!-- gretl customizations by Allin Cottrell -->

<style-sheet>

<style-specification id="print" use="docbook">
<style-specification-body> 

;; customize the print stylesheet

(define %generate-article-toc% 
  ;; Should a Table of Contents be produced for Articles?
  #t)

(define %generate-article-titlepage-on-separate-page%
  ;; Should the article title page be on a separate page?
  #t)

(define tex-backend 
  ;; Are we using the TeX backend?
  #t)

(define %default-quadding%   
  ;; The default quadding
  'justify)

(define bop-footnotes
  ;; Make "bottom-of-page" footnotes?
  #t)

;; Use #t for final version?
(define %footnote-ulinks%
  ;; Generate footnotes for ULinks?
  #f)

;; If set to #f, the table will span the entire page width.
(define %simplelist-column-width% 
  ;; Width of columns in tabular simple lists
  #f)

;; If #t, rules will be drawn before and after each Table.
(define %table-rules%
  ;; Specify rules before and after an Table
  #f)

;; If true, the URL of each ULink will appear in parenthesis after the 
;; text of the link. If the text of the link and the URL are identical, 
;; the parenthetical URL is suppressed.  
(define %show-ulinks%
  ;; Display URLs after ULinks?
  #f)

(define %hyphenation%
  ;; Allow automatic hyphenation?
  #t)

;;What font would you like for titles?
(define %title-font-family% 
  "LucidaSans")
  
(define %guilabel-font-family%
  "LucidaSans")
    
(define %admon-font-family%
  "LucidaSans")

;;What font would you like for the body?
(define %body-font-family% 
 "LucidaBright")

;;What font would you like for mono-seq?
(define %mono-font-family% 
 "LucidaTypewriter")

;;If the base fontsize is 10pt, and '%hsize-bump-factor%' is
;; 1.2, hsize 1 is 12pt, hsize 2 is 14.4pt, hsize 3 is 17.28pt, etc
(define %hsize-bump-factor% 
 1.1)
;;1.2

(define %ss-size-factor% 0.7)
(define %ss-shift-factor% 0.3)

(define (BULLSIZE m lvl)
  (let ((md (case-fold-down m)))
    (case md
          (("bullet") (MSIZE m lvl 1.0 0.8))
          (("box") (MSIZE m lvl 0.9 0.72))
          (("checkbox") (MSIZE m lvl 0.9 0.72))
          (("check") (MSIZE m lvl 1.0 1.0))
          (("checkedbox") (MSIZE m lvl 1.0 1.0))
          (("dash") (MSIZE m lvl 1.0 1.0))
          (("none") (MSIZE m lvl 1.0 1.0))
          (else (MSIZE m lvl 1.0 1.0)))))

(element (varlistentry term)
    (make paragraph
          space-before: (if (first-sibling?)
                            %block-sep%
                            0pt)
          keep-with-next?: #t
          first-line-start-indent: 0pt
          start-indent: 0pt
          font-posture: 'italic
          (process-children)))

;; These elements appear in this order on the title page of a book.
(define (book-titlepage-recto-elements)
  (list
        (normalize "title")
        (normalize "graphic")
        (normalize "subtitle")
        (normalize "author")
        (normalize "date")))
        
(define %footnote-size-factor% 
  ;; When printing footnotes, the current font size is multiplied by the
  ;; '%footnote-size-factor%'.  Default 0.9
  0.8)

(define %two-side% 
  ;; If '%two-side%' is true, headers and footers are alternated
  ;; so that the "outer" and "inner" headers will be correctly
  ;; placed in the bound document.  Default #f
  #t)

(define %admon-graphics%
  ;; If true, admonitions are presented in an alternate style that uses
  ;; a graphic.  Default graphics are provided in the distribution.
  #t)

(define %admon-graphics-path%
  ;; Sets the path, probably relative to the directory where the HTML
  ;; files are created, to the admonition graphics.
  "figures/")

(define admon-graphic-default-extension
  ;; Identifies the default extension for admonition graphics. This allows
  ;; backends to select different images (e.g., EPS for print, PNG for
  ;; PDF, etc.)
  ".png")

(define %line-spacing-factor% 
  ;; The leading is calculated by multiplying the current font size by the 
  ;; '%line-spacing-factor%'. For example, if the font size is 10pt and
  ;; the '%line-spacing-factor%' is 1.1, then the text will be
  ;; printed "10-on-11".  Default 1.3
  1.2)

(define %head-before-factor% 
  ;; The space before a title is calculated by multiplying the font size
  ;; used in the title by the '%head-before-factor%'.  Default 0.75
  0.4)

(define %head-after-factor% 
  ;; The space after a title is calculated by multiplying the font size used
  ;; in the title by the '%head-after-factor%'.  Default 0.5
  0.2)

(define %block-sep% 
  ;; The '%block-sep%' is the vertical distance between
  ;; block elements (figures, tables, etc.)
  ;; Default (* %para-sep% 2.0)
  (* %para-sep% 1.5))

(define formal-object-float
  ;; If '#t', formal objects will float if floating is supported by the
  ;; backend. At present, only the TeX backend supports floats.
  #f)

;;======================================
;;Inlines
;;======================================

(element application ($mono-seq$))
(element command ($mono-seq$))

&texmath.dsl;
&titlepage.dsl;

;; end of print stylesheet customization

</style-specification-body>
</style-specification>

<style-specification id="html" use="docbook">
<style-specification-body> 

;; customize the html stylesheet

;; This specifies the HTML extension to put on output files.
(define %html-ext% ".html")

;; Name for the root HTML document (default "book1")
(define %root-filename% "index")

(element application ($mono-seq$))
(element command ($mono-seq$))

&htmlmath.dsl;

;; end of html stylesheet customization

</style-specification-body>
</style-specification>

<external-specification id="docbook" document="docbook.dsl">

</style-sheet>

