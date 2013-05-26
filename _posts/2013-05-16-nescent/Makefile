index.html: header.html content.html footer.html
	cat $^ > index.html

content.html: content.md
	markdown $< > $@
