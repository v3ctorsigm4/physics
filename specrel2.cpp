

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <X11/Xlib.h>

typedef long double num;

num c = 3e8;

void lorentztran(num x1,num t1,num &x2,num &t2,num v)
{
	num gamma=1.0/sqrt(1.0-v*v/(c*c));

	x2=(x1-v*t1)*gamma;
	t2=(t1-v*x1/(c*c))*gamma;
}	


int main()
{
	num x1=0,t1=1.0/c,x2,t2,v,s,s2,s3,t3;
	
	v=0.75*c;
	s=t1*c;
	printf("%.12Le s\n",s);
	printf("%.12Le t1\n",t1);


	t2=t1*c/sqrt(c*c-v*v);
	printf("%.12Le t1b\n",t2);
	
//	s2=t2/(1/(c+v) + 1/(c*c/v -v));
//	s2=t2/(1/(c+v) + v/(c*c -v*v));
	s2=t2*2.0/(1/(c+v) + 1/(c -v));
	printf("%.12Le s1b\n",s2);

	
	t3=2/sqrt(c*c-v*v);
	printf("%.12Le t3\n",t3);
	

//	s3=s*sqrt(1.0-v*v/(c*c));
//	s3=t2*(c*c - v*v)/(2*c);//felut!
//	s3=2*sqrt(c*c-v*v)*(c*c - v*v)/(2*c);
	s3=2*(c*c - v*v)/(2*c*sqrt(c*c-v*v));
	printf("%.12Le s3\n",s3);

	
	t2=t1/sqrt(1.0-v*v/(c*c));
	printf("%.12Le t1a timedilation\n",t2);

	s3=s*sqrt(1.0-v*v/(c*c));
	printf("%.12Le s4 length contr\n",s3);





	lorentztran(x1,t1,x2,t2,v);
	printf("%.12Le t2 Lorentz\n",t2);

	s2=s*sqrt(1.0-v*v/(c*c));
	printf("%.12Le s2 Lorentz\n",s2);

	return 0;
}
/*

ORIGINAL OF SPECREL
hatra a fenyjel es a vegtelen jel metszete
t3 = t1+t2
t1 = s/(c+v)
t2 = s/(c*c/v -v) = t3-t1

t1=s/(c+v)
t1=t3-s/(c*c/v -v)
s/(c+v)=t3-s/(c*c/v -v)
s/(c+v) + s/(c*c/v -v)=t3
s*(1/(c+v) + 1/(c*c/v -v))=t3
s=t3/(1/(c+v) + 1/(c*c/v -v))

s=t3/(1/(c+v) + v/(c*c -v*v))
s=t3/(((c*c -v*v) + (c+v)v)/ ((c*c -v*v)(c+v)))
s=t3/(((c*c -v*v) + (c+v)v)/ ((c*c*c -v*v*c)+(c*c*v -v*v*v)))
s=t3/((c*c -v*v + cv+vv/ (c*c*c -v*v*c + c*c*v -v*v*v))
s=t3/((c*c+ cv/ (c*c*c -v*v*c + c*c*v -v*v*v))


!!!!!!!!!!!!!!!
elore-hatra a fenyjelek atlag erkezesi ideje az dilatalt ido!
t3=(t1+t2)/2

t1=s/(c+v) //hatra
t2=s/(c-v)= (2t3-t1)  //elore
(2t3-t1)=s/(c-v)
t1=2t3-s/(c-v)

s/(c+v)=2t3-s/(c-v)
s/(c+v)+s/(c-v)=2t3
s*(1/(c+v)+1/(c-v))=2t3
s=2t3/(1/(c+v)+1/(c-v))



?????????????????
elore
t=(s+t*v)/c
tc=s+t*v
tc-t*v=s
t*(c-v)=s
t=s/(c-v)

hatra
t=(s-t*v)/c
tc=s-t*v
tc+t*v=s
t*(c+v)=s
t=s/(c+v)


t1=s/(c+v)  
t2=s/(c-v)  

t3=t1+t2 t2=t3-t1
t3-t1=s/(c-v)  
t1=-(s/(c-v)-t3)
t1=t3-s/(c-v)

t1=s/(c+v)  
t1=t3-s/(c-v)
s/(c+v)=t3-s/(c-v)
s/(c+v)+s/(c-v)=t3
s*(1/(c+v) + 1/(c-v))=t3
s=t3/(1/(c+v) + 1/(c-v))
s=t3/(((c+v)+(c-v)) / ((c-v)*(c+v)))
s=t3/((c+c)/(c*c-v*c + c*v-v*v))
s=t3/(2c/(c*c - v*v))
s=t3*(c*c - v*v)/(2c)

s0=1
t0=2*s0/c=2/c   oda vissza ido
t3=t0*c/sqrt(c*c-v*v)
t3=2/sqrt(c*c-v*v)

s=2*(c*c - v*v)/(2c*sqrt(c*c-v*v))
s=sqrt(c*c - v*v)/c





  

Na, akkor kezdjük. Mi is a lényege a relativitásnak és honnan ered az idődilatáció és a hosszkontrakció vagyis a Minkowski geometria? Van egy rövid válasz erre, de csak és szigorúan a matek után.

Ugyanis az a fizika nyelve.


Szóval a fény sebessége független a forrás sebességétől. A következő levezetés erre az egyetlen tényre épít.

Legyen egy közegünk és mozogjon benne a hullám. 2/c másodperc a hullám periódusideje. A hullámhosszt úgy definiáljuk, hogy a hullám visszaverődik egy (x0) koordinátáról. A visszaverődő hullámok jól meghatározott téridő koordinátákban fogják metszeni egymást. Ezek a konstruktív interferencia helyek fogják megadni az adott "test" méretét.

Tehát egy álló IR (inercia rendszer)-ben a fentebb említett periódusidejű hullám hullámhossza 1 méter a faltól (x0) számítva. És innen periódikusan ismétlődik a konstruktív interferenciájú hely.

Tehát azért választottam 2/c periódusidőt, hogy a hullám vissza tudjon érni pont a következő hullámfront indulásáig, ha s=1 méter. Ez a hullámhossz az álló IR-ben.


Most mozogjon az egész cucc (v) sebességgel. Igy néz ki ez egy tér és egy idő dimenzióban. A koordináta rendszer ugyan az, mint amiben az előbb leírtam a folyamatot. Tehát nincs és nem is volt semmiféle koordináta transzformáció.

A piros fényjel előre és hátra az alábbi egyenletekkel számolható idő alatt ér el a kék mozgó x koordinátától a sárga x0 falig, ami most nyilván mozog.

előre mozgó fényjel futási ideje:
t=(s+t*v)/c
tc=s+t*v
tc-t*v=s
t*(c-v)=s
t=s/(c-v)

és hátra
t=(s-t*v)/c
tc=s-t*v
tc+t*v=s
t*(c+v)=s
t=s/(c+v)

Tehát a legalsó piros vonalnak az oda-vissza út a kék konstruktív interferencia helyig t3 időbe telik. Ahol

t1=s/(c+v)  
t2=s/(c-v)  
t3=t1+t2 t2=t3-t1


Ebből kiszámolható az (s) távolság függése a t3-tól és a (v) sebességtől. A (c) a hullám terjedési sebessége.

t3-t1=s/(c-v)  
t1=-(s/(c-v)-t3)
t1=t3-s/(c-v)

t1=s/(c+v)  
t1=t3-s/(c-v)
s/(c+v)=t3-s/(c-v)
s/(c+v)+s/(c-v)=t3
s*(1/(c+v) + 1/(c-v))=t3
s=t3/(1/(c+v) + 1/(c-v))
s=t3/(((c+v)+(c-v)) / ((c-v)*(c+v)))
s=t3/((c+c)/(c*c-v*c + c*v-v*v))
s=t3/(2c/(c*c - v*v))
s=t3*(c*c - v*v)/(2c)


Mostmár csak a t3 függése kellene. Erre Feynman kitalált egy jó modellt, a fényórát. A mozgásra merőlegesen pattogó fényjellel időt tudunk mérni. Ha ez a fényóra mozog, akkor a mért idő lassul. Az ok egyszerű: a jelnek keresztbe kell haladnia. Ebből a keresztmozgásból számolható az idődilatáció az alábbiak szerint

q=c/sqrt(c*c - v*v)  Ez a mozgásra merőleges sebesség komponens aránya a c-hez.

rendezzük át az egyenletet,
q=(c/c)/(sqrt(c*c - v*v) /c)
q=1/sqrt((c*c - v*v)/(c*c)) 
q=1/sqrt(1 - v*v/(c*c))
https://en.wikipedia.org/wiki/Time_dilation
15–4Transformation of time
http://www.feynmanlectures.caltech.edu/I_15.html
Az idődilatáció rendben.


Mostmár tudjuk a t3 értékét a mozgó IR-ben.

s0=1 méter
t0=2*s0/c=2/c   az oda-vissza idő az első IR-ben
t3=t0*c/sqrt(c*c-v*v) és a mozgó IR-ben
t3=2/sqrt(c*c-v*v)

s=2*(c*c - v*v)/(2c*sqrt(c*c-v*v))
s=sqrt(c*c - v*v)/c

Ez a hosszkontrakció egyenlete, ha a test méretét egy hullám konstruktív interferencia helyei határozzák meg. Rendezzük át az egyenletet.

s=sqrt((c*c - v*v)/(c*c))
s=sqrt(1 - v*v/(c*c))
s0=1
s=1sqrt(1 - v*v/(c*c))
https://en.wikipedia.org/wiki/Length_contraction

Szóval a rövid válasz:
azért ilyen a "geometriája" a világnak, mert minden csak közönséges hullám.

 

 







Coherent Excited States in the Theory of Superconductivity:
Invariance and the Meissner Effect


https://en.wikipedia.org/wiki/Weinberg_angle
Tehát az elektromágneses tér két komponensű, egy (B) skalármező és egy (W) vektormező keveréke.
De hát ezt eddig is így tudtuk, csak más volt a neve a gyereknek.
https://en.wikipedia.org/wiki/Electromagnetic_four-potential
Skalár-potenciál és vektor-potenciál.

 

A (B) skalármező az nem a Higgs mező, hanem a hypercharge-mező. Igazából a kettő szoros kapcsolatban áll, hiszen a Higgs-mező várhatóértéke nemzéró a vákuumban. Magyarul a vákumnak hypercharge töltése van. Ez olyasmi, mint a Dirac-tenger de mégsem az, hiszen ez a töltés nem elektromos töltés. A Dirac-tenger a félvezetők helyes modellje.

 

A foton-mező transzverzális rezgés, amit kétféle képpen lehet "kikeverni" ezekből a terekből. Az egyik módja az, hogy csak a vektorpoteciált forgatom a menetirányra merőlegesen. Ekkor ennek a vektornak a deriváltja kiadja a elektromos mező vektorát. Ez a E = -dA / dt tag. A mágneses mező ennek a curlja, B=nabla X A.

Ellenben a foton-tér a két összetevő keveréke, mint láttuk. Tehát a következő megoldás lesz a helyes.

A vektorpotenciál nem merőlegesen forog a menetirányra, hanem csak egy komponense merőleges arra. A másik komponense a menetirányba esik. DE a skalármező gradiense pont ellentétes irányú lesz a -dA / dt tag menetirányba eső komponensével, ezért az eltűnik.

A következő három ábra az utóbbi esettel számoltam. Az első ábrán a kék az (E) elektromos tér vektora a piros a (B) mágneses tér vektora. Merőlegesek egymásra és a mozgásirányra, ahogy illik.

A második ábrán a kék vektorok a skalármző gradiens vektorai. A harmadik ábra a "spinwave" avagy  maga az (A) vektormező.

*/


