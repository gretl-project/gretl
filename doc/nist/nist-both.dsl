<!DOCTYPE style-sheet PUBLIC "-//James Clark//DTD DSSSL Style Sheet//EN" [
<!ENTITY % html "IGNORE">
<![%html;[
<!ENTITY % print "IGNORE">
<!ENTITY docbook.dsl SYSTEM "/usr/share/sgml/docbook-dsssl-1.74b/html/docbook.dsl" CDATA dsssl>
]]>
<!ENTITY % print "INCLUDE">
<![%print;[
<!ENTITY docbook.dsl SYSTEM "/usr/share/sgml/docbook-dsssl-1.74b/print/docbook.dsl" CDATA dsssl>
]]>
<!ENTITY htmlmath.dsl SYSTEM "../HTMLMath.dsl">
<!ENTITY texmath.dsl SYSTEM "../TeXMath.dsl">
]>

<!-- customizations by Allin Cottrell -->

<style-sheet>

<style-specification id="print" use="docbook">
<style-specification-body> 

;; customize the print stylesheet

(define %generate-article-toc% 
  ;; Should a Table of Contents be produced for Articles?
  #f)

(define %generate-article-titlepage-on-separate-page%
  ;; Should the article title page be on a separate page?
  #f)

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
          font-family-name: %mono-font-family%
          (process-children)))

;; These elements appear in this order on the title page of a book.
(define (article-titlepage-recto-elements)
  (list
        (normalize "title")
        (normalize "author")
        (normalize "pubdate")))

(define %footnote-size-factor% 
  ;; When printing footnotes, the current font size is multiplied by the
  ;; '%footnote-size-factor%'.  Default 0.9
  0.8)

(define %line-spacing-factor% 
  ;; The leading is calculated by multiplying the current font size by the 
  ;; '%line-spacing-factor%'. For example, if the font size is 10pt and
  ;; the '%line-spacing-factor%' is 1.1, then the text will be
  ;; printed "10-on-11".  Default 1.3
  1.2)

(define %head-before-factor% 
  ;; The space before a title is calculated by multiplying the font size
  ;; used in the title by the '%head-before-factor%'.  Default 0.75
  0.5)

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

;; end of print stylesheet customization

</style-specification-body>
</style-specification>

<style-specification id="html" use="docbook">
<style-specification-body> 

;; customize the html stylesheet

;; This specifies the HTML extension to put on output files.
(define %html-ext% ".html")

(define %generate-article-toc% 
  ;; Should a Table of Contents be produced for Articles?
  ;; If true, a Table of Contents will be generated for each 'Article'.
  #t)

;; put keywords onto HTML title page
(define (article-titlepage-recto-elements)
  (list
        (normalize "title")
        (normalize "author")
        (normalize "pubdate")
        (normalize "keywordset")))

(mode article-titlepage-recto-mode

  (element affiliation
    (make element gi: "DIV"
	  attributes: (list (list "CLASS" (gi)))
	  (process-children)))

  (element author
    (let ((author-name  (author-string))
	  (author-affil (select-elements (children (current-node)) 
					 (normalize "affiliation"))))
      (make sequence      
	(make element gi: "H3"
	      attributes: (list (list "CLASS" (gi)))
	      (make element gi: "A"
		    attributes: (list (list "NAME" (element-id)))
		    (literal author-name)))
	(process-node-list author-affil))))

  (element date
    (make element gi: "SPAN"
	  attributes: (list (list "CLASS" (gi)))
	  (process-children)
	  (make empty-element gi: "BR")))

  (element firstname
    (make element gi: "SPAN"
	  attributes: (list (list "CLASS" (gi)))
	  (process-children)
	  (make empty-element gi: "BR")))

  (element graphic
    (let* ((nd (current-node))
	   (fileref (attribute-string (normalize "fileref") nd))
	   (entattr (attribute-string (normalize "entityref") nd))
	   (entityref (if entattr
			  (entity-system-id entattr)
			  #f))
	   (format  (attribute-string (normalize "format")))
	   (align   (attribute-string (normalize "align")))
	   (attr    (append 
		     (if align 
			 (list (list "ALIGN" align)) 
			 '())
		     (if entityref
			 (list (list "SRC" (graphic-file entityref)))
			 (list (list "SRC" (graphic-file fileref))))
		     (list (list "ALT" ""))
		     )))
      (if (or fileref entityref) 
	  (make empty-element gi: "IMG"
		attributes: attr)
	  (empty-sosofo))))

  (element keywordset 
       (make element gi: "P" 
           attributes: (list (list "CLASS" (gi)))
           (literal "Keywords: ")
           (process-children)
           (make empty-element gi: "BR")))
           
  (element keyword
     (if (last-sibling?) 
        (literal (data (current-node))) 
        (literal (string-append (data (current-node)) ", "))))
	
  (element legalnotice 
    (titlepage-recto-legalnotice))
  
  (element (legalnotice title) (empty-sosofo))

  (element orgdiv
    (make element gi: "SPAN"
	  attributes: (list (list "CLASS" (gi)))
	  (process-children)
	  (make empty-element gi: "BR")))

  (element orgname
    (make element gi: "SPAN"
	  attributes: (list (list "CLASS" (gi)))
	  (process-children)
	  (make empty-element gi: "BR")))

  (element pubdate
    (make element gi: "P"
	  attributes: (list (list "CLASS" (gi)))
	  (process-children)
	  (make empty-element gi: "BR")))

  (element releaseinfo
    (make element gi: "SPAN"
	  attributes: (list (list "CLASS" (gi)))
	  (process-children)
	  (make empty-element gi: "BR")))

  (element subtitle 
    (make element gi: "H2"
	  attributes: (list (list "CLASS" (gi)))
	  (process-children-trim)))

  (element surname
    (make element gi: "SPAN"
	  attributes: (list (list "CLASS" (gi)))
	  (process-children)
	  (make empty-element gi: "BR")))

  (element title 
    (make element gi: "H1"
	  attributes: (list (list "CLASS" (gi)))
	  (make element gi: "A"
		attributes: (list (list "NAME" (element-id)))
		(with-mode title-mode
		  (process-children-trim)))))

  (element titleabbrev
    (make element gi: "SPAN"
	  attributes: (list (list "CLASS" (gi)))
	  (process-children)
	  (make empty-element gi: "BR")))
)

(define %shade-verbatim% #t)

(define %section-autolabel% #t)

(define ($section-separator$)
  (empty-sosofo)) 

(define %stylesheet% "dbtex.css")

(define %html-pubid%
  ;; REFENTRY html-pubid
  ;; PURP What public ID are you declaring your HTML compliant with?
  ;; DESC
  ;; The public ID used in output HTML files.  If '#f', then no doctype
  ;; declaration is produced.
  ;; /DESC
  ;; AUTHOR N/A
  ;; /REFENTRY
  "-//W3C//DTD HTML 4.0 Transitional//EN")

(define %html-header-tags% 
  ;; REFENTRY html-header-tags
  ;; PURP What additional HEAD tags should be generated?
  ;; DESC
  ;; A list of the the HTML HEAD tags that should be generated.
  ;; The format is a list of lists, each interior list consists
  ;; of a tag name and a set of attribute/value pairs:
  ;; '(("META" ("NAME" "name") ("CONTENT" "content")))
  ;; /DESC
  ;; AUTHOR N/A
  ;; /REFENTRY
 '(("META" ("http-equiv" "Content-Type") ("CONTENT" "text/html; charset=utf-8"))))

(define ($graphic$ fileref
                   #!optional (format #f) (alt #f) (align #f) 
		   (width #f) (height #f))
  (let ((img-attr  (append
                    (list     (list "SRC" (graphic-file fileref)))
                    (if alt   (list (list "ALT" alt)) '())
                    (if align
			(if (equal? align "CENTER") 
			(list (list "ALIGN" "MIDDLE"))
			(list (list "ALIGN" align)))'())
                    (if width (list (list "WIDTH" width)) '())
                    (if height (list (list "HEIGHT" height)) '()))))
    (make empty-element gi: "IMG"
          attributes: img-attr)))

(element application ($mono-seq$))
(element command ($mono-seq$))

&htmlmath.dsl;

;; end of html stylesheet customization

</style-specification-body>
</style-specification>

<external-specification id="docbook" document="docbook.dsl">

</style-sheet>

