TFYA50
=====

Gör:
-När du börjar för dagen är det en bra vana att hämta senaste versionen av alla
filer: 
   $ git pull

-När du är klar med något OCH KODEN FUNGERAR(!!!) så är det dags att committa
 ändringarna, dvs ta ett "foto" av hur filerna ser ut:
   $ git add <filnamn>
   ... (Gör för alla filer du ändrat i)
   $ git commit -m "Beskriv vad du gjort"
   
-Nu när du har tagit ditt foto så ska du ladda upp det till servern (push)

   $ git push

(Det kan hända att git klagar på att du inte har senaste versionen av projektet.
 I så fall så behöver du köra 'pull' en gång till innan du pushar)

Gör INTE:
- Gör inte commit eller push med filer som inte fungerar. Snälla ladda inte upp
 filer till servern om de orsakar krachar eller får programmet att uppföra sig
 konstigt. Vänta med att ladda upp dem tills du är klar med det du håller på
 med. 

