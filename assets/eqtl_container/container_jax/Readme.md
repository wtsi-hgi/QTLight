Note that we have 
git clone https://github.com/mancusolab/jaxqtl.git
commit c126cbc381417eab7a33fd16cf4635adfde4ba1c (HEAD -> main, origin/main, origin/dev, origin/HEAD)

And then we patched the JaxQTL code to allow numeric sample identifiers:
have to make sure the covar index is string. I added covar.index = covar.index.astype(str)  to jaxqtl/src/jaxqtl/io/readfile.py  line 176 so before taking the intersection. you can avoid that by using sample ids which cannot be confused with integers