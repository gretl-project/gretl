<!DOCTYPE style-sheet PUBLIC "-//James Clark//DTD DSSSL Style Sheet//EN" [
<!ENTITY dbstyle SYSTEM "/usr/share/sgml/docbook-dsssl-1.72/print/docbook.dsl" 
  CDATA DSSSL>
<!ENTITY mathml.dsl SYSTEM "MathMLrev.dsl">
<!ENTITY titlepage.dsl SYSTEM "titlepage.dsl">
]>

<style-sheet>
<style-specification use="docbook">
<style-specification-body>

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
  #t)

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
  "Helvetica")

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

;; These elements appear in this order on the title page of a book.
(define (book-titlepage-recto-elements)
  (list
        (normalize "title")
	(normalize "graphic")
        (normalize "subtitle")
        (normalize "author")
        (normalize "date")))


;;======================================
;;Inlines
;;======================================

(element application ($mono-seq$))
(element command ($mono-seq$))

&mathml.dsl;
&titlepage.dsl;

</style-specification-body>
</style-specification>
<external-specification id="docbook" document="dbstyle">
</style-sheet>