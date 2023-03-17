#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <sys/stat.h>
#include <math.h> 
// ceil, floor : #include <math.h>
#include "./include/traces.h" 
#include "./include/check.h" 
 

#define isINTREE(i,n)		(i<n)
#define iPARENT(i) 			(i-1)/2
#define iLCHILD(i) 			(2*i)+1
#define iRCHILD(i) 			(2*i)+2



#define MAXCARS 128

typedef char T_elt;


typedef struct {
	unsigned int nbElt;
	unsigned char tree[MAXCARS]; //tableau des caractères
	int data[2*MAXCARS-1]; //tableau des occurences
	int huffmanTree[2*MAXCARS -1];
} T_indirectHeap;


//Macro-fonctions
#define iPARENT(i) 			(i-1)/2
#define iLCHILD(i) 			(2*i)+1
#define iRCHILD(i) 			(2*i)+2
#define isLEAF(i,n) 			(2*i>=(n-1))
#define isINTREE(i,n)		(i<n)
#define isROOT(i)				(i==0)


#define VALP(pHeap, i)		pHeap->tree[i]		
#define VAL(heap, i)		heap.tree[i]
#define OCC(heap, i)		heap.data[i]
#define OCCP(heap, i)		heap->data[i]

//Prototypes
static void genDotPOT_rec(T_elt t[], int n, int root, FILE *fp); 
void createDotPOT(T_elt t [], int n, const char *basename); 
void createDotPOT_rec(T_elt t [], int n, int root, const char *basename, int indent);
T_indirectHeap * newHeap();
T_indirectHeap *analyserDocument(char *document);
void printHeap(T_indirectHeap *);
void swap(T_indirectHeap *p, int i, int j);
void siftDown(T_indirectHeap *p, int k);
void buildHeapV2(T_indirectHeap *p);//transformerEnMinimierV2
void initHuffmanTree(int* huffmanTree);
T_elt removeMin(T_indirectHeap *p);
void insererMI(T_indirectHeap *p, int e, int occ);
void addNode(int* , int , int );
void huffman(char *document);
void printHuffmanTree(T_indirectHeap *p);




char * outputPath = "./tp5";

int main(void) {
	/*
	char * document = "ABRADACABRA";
	T_indirectHeap *p = analyserDocument(document);
	printHeap(p);
	buildHeapV2(p);
	printHeap(p);
	*/
	huffman("ABRADACABRA");
	//createDotPOT(p->tree, p->nbElt, "tas");
	return 0;
}

int eltCmp(T_elt e1, T_elt e2) {
	return e1-e2;
}


T_indirectHeap * newHeap(){
	T_indirectHeap * pAux= (T_indirectHeap *) malloc(sizeof(T_indirectHeap));
	pAux->nbElt = 0;
	for (int i = 0; i < 2*MAXCARS -1; i++) {
		pAux->huffmanTree[i] = -256;
	}
	return pAux;
}




static void genDotPOT_rec(T_elt t[], int n, int root, FILE *fp){
	// Attention : les fonction toString utilisent un buffer alloué comme une variable statique 
	// => elles renvoient toujours la même adresse 
	// => on ne peut pas faire deux appels à toString dans le même printf()
	
	// t : tas
	// n : taille du tas
	// root : indice de la racine du sous-arbre à produire
	// fp : flux correspondant à un fichier ouvert en écriture où écrire le sous-arbre

	// ordre de récurrence  ?

	// cas de base 
	if(!isINTREE(root,n)) return;

	// cas général 
	fprintf(fp, "\t%d",root); 
    fprintf(fp, " [label = \"%d\"];\n",(t[root]));
    if ( !isINTREE(iRCHILD(root),n) && !isINTREE(iLCHILD(root),n)) {
        fprintf(fp, "\t%d", root); 
		fprintf(fp, " [label = \"%c\"];\n", (t[root]));
	}
    else if (!isINTREE(iRCHILD(root),n)) {
        fprintf(fp, "\t%d", root);
		fprintf(fp, " [label = \"%d\"];\n", (t[root]));
	}
	else if (!isINTREE(iLCHILD(root),n)) {
        fprintf(fp, "\t%d",root);
		fprintf(fp, " [label = \"%d\"];\n", (t[root]));
	}
	
    if (isINTREE(iLCHILD(root),n)) {
        fprintf(fp, "\t%d",root);
		fprintf(fp, ":sw -> %d;\n", iLCHILD(root));
		genDotPOT_rec(t,n,iLCHILD(root),fp);
    }

    if (isINTREE(iRCHILD(root),n)) {
        fprintf(fp, "\t%d",root);
		fprintf(fp,":se -> %d;\n", iRCHILD(root));
		genDotPOT_rec(t,n,iRCHILD(root),fp);
    }

}


void createDotPOT(T_elt t [], int n,const char *basename) {
	static char oldBasename[FILENAME_MAX + 1] = "";
	static unsigned int noVersion = 0;

	char DOSSIER_DOT[FILENAME_MAX + 1]; 
	char DOSSIER_PNG[FILENAME_MAX + 1]; 

	char fnameDot [FILENAME_MAX + 1];
	char fnamePng [FILENAME_MAX + 1];
	char	cmdLine [2 * FILENAME_MAX + 20];
	FILE *fp;
	struct stat sb;
	

	// Au premier appel, création (si nécessaire) des répertoires
	// où seront rangés les fichiers .dot et .png générés par cette fonction	

	// il faut créer le répertoire outputPath s'il n'existe pas 
	if (stat(outputPath, &sb) == 0 && S_ISDIR(sb.st_mode)) {
    } else {
        printf("Création du répertoire %s\n", outputPath);
		mkdir(outputPath, 0777);
    }

	// il faut créer les répertoires outputPath/png et /dot 
	sprintf(DOSSIER_DOT, "%s/dot/",outputPath);
	sprintf(DOSSIER_PNG, "%s/png/",outputPath);

	if (oldBasename[0] == '\0') {
		mkdir(DOSSIER_DOT,	S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
		mkdir(DOSSIER_PNG,	S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
	}

	 // S'il y a changement de nom de base alors recommencer à zéro
	 // la numérotation des fichiers 

	if (strcmp(oldBasename, basename) != 0) {
		noVersion = 0;
		strcpy(oldBasename, basename); 
	}

	sprintf(fnameDot, "%s%s_v%02u.dot", DOSSIER_DOT, basename, noVersion);
	sprintf(fnamePng, "%s%s_v%02u.png", DOSSIER_PNG, basename, noVersion);

	CHECK_IF(fp = fopen(fnameDot, "w"), NULL, "erreur fopen dans saveDotBST"); 
	
	noVersion ++;
	fprintf(fp, "digraph %s {\n", basename);
	fprintf(fp, 
	"\tnode [\n"
		"\t\tfontname  = \"Arial bold\" \n"
		"\t\tfontsize  = \"14\"\n"
		"\t\tfontcolor = \"red\"\n"
		"\t\tstyle     = \"rounded, filled\"\n"
		"\t\tshape     = \"circle\"\n"
		"\t\tfillcolor = \"grey90\"\n"
		"\t\tcolor     = \"blue\"\n"
		"\t\twidth     = \"0.5\"\n"
		"\t]\n"
	"\n"
	"\tedge [\n"
		"\t\tcolor     = \"blue\"\n"
	"\t]\n\n"
	);


	genDotPOT_rec(t,n,0,fp);

	fprintf(fp, "}\n");	
	fclose(fp);

	sprintf(cmdLine, "dot -Tpng  %s -o %s", fnameDot, fnamePng);
	system(cmdLine);

	printf("Creation de '%s' et '%s' ... effectuee\n", fnameDot, fnamePng);
}


/**
* Mi : minimier indirect
* Ht : arbre de codage de Huffman
* C1, C2 : caractères extraits
* Ni : noeud inséré


huffman(D)
  Mi = analyserDocument(D)	// Comptage des occurrences
  Ht = initHuffmanTree() // Initialisation de l’arbre de codage
  transformerEnMinimierV2(Mi) // Réorganisation du minimier
  Pour i = 1 jusque n-1 			
	C1=extraireMin(Mi) // Extraire et réorganiser
	C2=extraireMin(Mi) // Extraire et réorganiser
	AjouterNoeud(Ht,Ni) // Ajout dans l’arbre de codage
	insererMI(Mi,Ni,VAL(C1)+VAL(C2)) // Insérer et réorganiser

*/

void huffman(char *document) {
	T_indirectHeap *Mi = analyserDocument(document);
	int C1, C2, n;
	buildHeapV2(Mi);

	printf("Creation de l'arbre...\n");
	printf("nbElt = %d\n", Mi->nbElt);
	printHeap(Mi);
	printf("\n");
	n= Mi->nbElt;
	for (int i = 1; i < n-1; i++) {
		C1 = removeMin(Mi);
		C2 = removeMin(Mi);
		Mi->huffmanTree[C2] = -(MAXCARS+i);
		Mi->huffmanTree[C1] = (MAXCARS+i);
		insererMI(Mi,(MAXCARS+i), OCCP(Mi,C1) + OCCP(Mi,C2));
	}
	printf("Arbre de codage de Huffman :\n");
	printHuffmanTree(Mi);
}



T_indirectHeap *analyserDocument(char *document) {
	T_indirectHeap *Mi = newHeap();
	int trouve = 0, intC;
	for (int i = 0; i < strlen(document) ; i++) {
		trouve = 0,intC=document[i];
		// on parcourt le tableau des caractères pour incrémenter le compteur et s'il est présent
		for (int j=0; j < Mi->nbElt; j++) {
			if (document[i] == Mi->tree[j]) {
				trouve = 1;
				break;
			}
		}
		if (trouve==0) {
			Mi->tree[Mi->nbElt] = document[i];
			Mi->nbElt++;
		}
		OCCP(Mi,intC) ++;
		if (Mi->tree[Mi->nbElt -1] == 9)
			printf("\n%d\n", Mi->tree[Mi->nbElt -1]);
		
	}
	return Mi;
}


void printHeap(T_indirectHeap *Mi) {
	printf("Nombre d'éléments : %d\n", Mi->nbElt);
	for (int i = 0; i < Mi->nbElt; i++) {
		printf("%d : %d\n", Mi->tree[i], OCCP(Mi, Mi->tree[i]) );
	}
}

void swap(T_indirectHeap *p, int i, int j) {
	T_elt aux = VALP(p,i);
	VALP(p,i) = VALP(p,j);
	VALP(p,j) = aux;
}


void siftDown(T_indirectHeap *p, int k) {
	int i, fini =  isLEAF(k,p->nbElt);
	while (!fini) {
		i = iLCHILD(k);
		if (isINTREE(iRCHILD(k), p->nbElt) && OCCP(p,VALP(p,iRCHILD(k))) < OCCP(p,VALP(p,iLCHILD(k))))
			i = iRCHILD(k);
		if (OCCP(p,VALP(p,k)) < OCCP(p,VALP(p,i))) 
			fini = 1;
		else {
			swap(p, k, i);
			k = i;
			fini = isLEAF(k,p->nbElt);
		}
	}
}

void buildHeapV2(T_indirectHeap * p){
	for (int i = p->nbElt/2; i >= 0; i--)
		siftDown(p, i);
}

T_elt removeMin(T_indirectHeap *p) {
	swap(p, 0, p->nbElt-1);
	p->nbElt--;
	siftDown(p, 0);
	return VALP(p,p->nbElt); 
}


void insererMI(T_indirectHeap *p, int e, int occ){
	VALP(p,p->nbElt) = e;
	p->nbElt++;
	OCCP(p,e) = occ;
	buildHeapV2(p);
}

void printHuffmanTree(T_indirectHeap *p) {
	int fini = 0;
	while (!fini)
	{
		for (int i = 0; i < p->nbElt; i++) {
			if (p->huffmanTree[i] < 0)
				printf("%d : %d\n", p->tree[i], p->huffmanTree[p->tree[i]]);
		}
		fini = 1;
	}
	
}