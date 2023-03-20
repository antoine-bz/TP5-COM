#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <sys/stat.h>
#include <math.h> 
// ceil, floor : #include <math.h>
#include "./include/traces.h" 
#include "./include/check.h" 

#define MAXCARS 128
#define ILASTCAR MAXCARS-1
#define IMAXCODAGE 16

typedef unsigned char T_elt;

typedef struct {
	unsigned int nbElt;
	T_elt tree[MAXCARS]; //tableau des caractères
	int data[2*MAXCARS-1]; //tableau des occurences
	int huffmanTree[2*MAXCARS -1];
	//on veut maintenant stocker
	char codage[MAXCARS][IMAXCODAGE];
} T_indirectHeap;


//Macro-fonctions
#define iPARENT(i) 			(i-1)/2
#define iLCHILD(i) 			(2*i)+1
#define iRCHILD(i) 			(2*i)+2
#define isLEAF(i,n) 		(2*i>=(n-1))
#define isINTREE(i,n)		(i<n)
#define isROOT(i)			(i==0)
#define nbINTERNALS(n) 		n/2


#define TREEP(pHeap, i)		pHeap->tree[i]		
#define TREE(heap, i)		heap.tree[i]
#define OCC(heap, i)		heap.data[i]
#define VALP(heap, i)		heap->data[i]

//Prototypes
static void genDotPOThuff_rec(T_indirectHeap *p, int root, FILE *fp); 
static void genDotPOTtree_rec(T_indirectHeap *p, int root, FILE *fp);
void createDotPOT(T_indirectHeap *p,const char *basename, int huffmanTree); 
T_indirectHeap * newHeap();
T_indirectHeap *analyserDocument(char *document);
void swapTree(T_indirectHeap *p, int i, int j);
void descendre(T_indirectHeap *p, int k);
void buildHeapV2(T_indirectHeap *p);//transformerEnMinimierV2
void insererMI(T_indirectHeap *p, int e, int occ);
void huffman(char *document);
void heapSortV2(T_indirectHeap *p);
void freeHeap(T_indirectHeap *p);
T_elt extraireMin(T_indirectHeap *p);

void printHeap(T_indirectHeap *);
void printHuffmanTree(T_indirectHeap *p);

void codageHuffman(T_indirectHeap *p);
void printCodage(T_indirectHeap *p);


char * outputPath = "./tp5";

int main(void) {
	/*
	char * document = "ABRADACABRA";
	T_indirectHeap *p = analyserDocument(document);
	printHeap(p);
	buildHeapV2(p);
	printHeap(p);
	*/
	//huffman("ABRADACABRA");
	huffman("algorithme de huffman pour la compression de chaines");
	//createDotPOT(p->tree, p->nbElt, "tas");
	return 0;
}


/**
 * nom : analyserDocument
 * but : analyser un document et retourner un tas contenant les caractères et leur 
 * @return un tas contenant les caractères et leur occurence
 */
T_indirectHeap * newHeap(){
	T_indirectHeap * pAux= (T_indirectHeap *) malloc(sizeof(T_indirectHeap));
	pAux->nbElt = 0;
	for (int i = 0; i < 2*MAXCARS -1; i++) {
		pAux->huffmanTree[i] = -256;
	}
	return pAux;
}

void freeHeap(T_indirectHeap *p){
	free(p);
}


static void genDotPOThuff_rec(T_indirectHeap *p, int root, FILE *fp){
	// Attention : les fonction toString utilisent un buffer alloué comme une variable statique 
	// => elles renvoient toujours la même adresse 
	// => on ne peut pas faire deux appels à toString dans le même printf()
	
	// t : tas
	// n : taille du tas
	// root : indice de la racine du sous-arbre à produire
	// fp : flux correspondant à un fichier ouvert en écriture où écrire le sous-arbre

	int i = root, i1=-256, i2=-256;


	for(int j = 0; j < MAXCARS*2; j++) {
		if(p->huffmanTree[j] == -i) {
			fprintf(fp, "\t%d", root);
			fprintf(fp, " [label = \"%d\"];\n", (root));
			fprintf(fp, "\t%d",root);
			fprintf(fp, ":sw -> %d [label= \"0\"];\n", j);
			if (j <= ILASTCAR)
				fprintf(fp, "\t%d [label = \"%c\"];\n",j,j);
			i1 = j;
		}
		if(p->huffmanTree[j] == i) {
			fprintf(fp, "\t%d", root);
			fprintf(fp, " [label = \"%d\"];\n", (root));
			fprintf(fp, "\t%d",root);
			fprintf(fp,":se -> %d [label= \"1\"];\n", j);
			if (j <= ILASTCAR)
				fprintf(fp, "\t%d [label = \"%c\"];\n",j,j);
			i2=j;
		}
	}
	if(i1 != i2) {
		if(i1 > ILASTCAR)
			genDotPOThuff_rec(p, i1, fp);
		if(i2 > ILASTCAR)
			genDotPOThuff_rec(p, i2, fp);
	}
	else 
		if(i1 >= ILASTCAR)
			genDotPOThuff_rec(p, i1, fp);

}

static void genDotPOTtree_rec(T_indirectHeap *p, int root, FILE *fp){
	// Attention : les fonction toString utilisent un buffer alloué comme une variable statique 
	// => elles renvoient toujours la même adresse 
	// => on ne peut pas faire deux appels à toString dans le même printf()
	
	// t : tas
	// n : taille du tas
	// root : indice de la racine du sous-arbre à produire
	// fp : flux correspondant à un fichier ouvert en écriture où écrire le sous-arbre

	// ordre de récurrence  ?

	//printHeap(p);

	int n = p->nbElt;

	// cas de base 
	if(!isINTREE(root,n))return;

	// cas général 
	fprintf(fp, "\t%d",root); 
    fprintf(fp, " [label = <%c<BR/><FONT POINT-SIZE=\"10\">%d</FONT>>];\n",TREEP(p,root),VALP(p,TREEP(p,root)));
    if ( !isINTREE(iRCHILD(root),n) && !isINTREE(iLCHILD(root),n)) {
        fprintf(fp, "\t%d", root); 
		fprintf(fp, " [label = <%c<BR/><FONT POINT-SIZE=\"10\">%d</FONT>>];\n",TREEP(p,root),VALP(p,TREEP(p,root)));
	}
    else if (!isINTREE(iRCHILD(root),n)) {
        fprintf(fp, "\t%d", root);
		fprintf(fp, " [label = <%c<BR/><FONT POINT-SIZE=\"10\">%d</FONT>>];\n",TREEP(p,root),VALP(p,TREEP(p,root)));
	}
	else if (!isINTREE(iLCHILD(root),n)) {
        fprintf(fp, "\t%d",root);
		fprintf(fp, " [label = <%c<BR/><FONT POINT-SIZE=\"10\">%d</FONT>>];\n",TREEP(p,root),VALP(p,TREEP(p,root)));
	}
	
    if (isINTREE(iLCHILD(root),n)) {
        fprintf(fp, "\t%d",root);
		fprintf(fp, ":sw -> %d [arrowhead = \"none\"];\n", iLCHILD(root));
		genDotPOTtree_rec(p,iLCHILD(root),fp);
    }

    if (isINTREE(iRCHILD(root),n)) {
        fprintf(fp, "\t%d",root);
		fprintf(fp,":se -> %d [arrowhead = \"none\"];\n", iRCHILD(root));
		genDotPOTtree_rec(p,iRCHILD(root),fp);
    }

}


void createDotPOT(T_indirectHeap *p,const char *basename, int huffmanTree) {
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
		"\t\tfontcolor = \"black\"\n"
		"\t\tstyle     = \"rounded, filled\"\n"
		"\t\tshape     = \"circle\"\n"
		"\t\tfillcolor = \"grey90\"\n"
		"\t\tcolor     = \"black\"\n"
		"\t\twidth     = \"0.5\"\n"
		"\t]\n"
	"\n"
	"\tedge [\n"
		"\t\tcolor     = \"black\"\n"
	"\t]\n\n"
	);


	if (huffmanTree == 1)
		genDotPOThuff_rec(p,p->tree[0],fp);
	else{
		genDotPOTtree_rec(p,0,fp);
	}

	fprintf(fp, "}\n");	
	fclose(fp);

	sprintf(cmdLine, "dot -Tpng  %s -o %s", fnameDot, fnamePng);
	system(cmdLine);

	printf("Creation de '%s' et '%s' ... effectuee\n", fnameDot, fnamePng);
}

/**
 * nom : huffman
 * but : créer l'arbre de codage de Huffman
 * @param document : nom du fichier à analyse
 * @return : void
*/
void huffman(char *document) {
	int C1, C2, n, Ni;

	T_indirectHeap *Mi = analyserDocument(document);
	heapSortV2(Mi);

	printf("Creation de l'arbre...\n");
	// printf("nbElt = %d\n", Mi->nbElt);
	// printHeap(Mi);
	printf("\n");

	createDotPOT(Mi, "tas",0);

	n= Mi->nbElt;
	for (int i = 1; i <= n-1; i++) {

		//extraireTop
		C1 = extraireMin(Mi);

		//extraireTop
		C2 = extraireMin(Mi);

		//AjouterNoeud
		Ni=ILASTCAR+i;
		Mi->huffmanTree[C1] = -Ni;
		Mi->huffmanTree[C2] = Ni;

		//insererMI
		insererMI(Mi,Ni, VALP(Mi,C1) + VALP(Mi,C2));

		// printf("\naffichage de l'arbre après l'insertion du noeud %d\n", Ni);
		// printHuffmanTree(Mi);
		// printf("\n\n");
	}
	createDotPOT(Mi, "huffmanTree",1);
	codageHuffman(Mi);
	printCodage(Mi);
	freeHeap(Mi);
}


/**
 * nom : analyserDocument
 * but : compter le nombre d'occurences de chaque caractère dans le document et les stocker les données dans les tableaux
 * @param document : nom du fichier à analyse
 * @return : T_indirectHeap *Mi
*/
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
		VALP(Mi,intC) ++;
		if (Mi->tree[Mi->nbElt -1] == 9)
			printf("\n%d\n", TREEP(Mi,Mi->nbElt -1));
		
	}
	return Mi;
}

/**
 * nom : printHeap
 * but : afficher le contenu du minimier
 * @param Mi : le minimier
 * @return void
*/
void printHeap(T_indirectHeap *Mi) {
	printf("Nombre d'éléments : %d\n", Mi->nbElt);
	for (int i = 0; i < Mi->nbElt; i++) {
		printf("%d : %d\n", Mi->tree[i], VALP(Mi, Mi->tree[i]) );
	}
}

/**
 * nom : swapTree
 * but : échanger deux éléments dans l'arbre
 * @param p : l'arbre
 * @param i : l'indice de l'élément à échanger
 * @param j : l'indice de l'élément à échanger
 * @return void
*/
void swapTree(T_indirectHeap *p, int i, int j) {
	if(i==j) return;
	T_elt aux = TREEP(p,i);
	TREEP(p,i) = TREEP(p,j);
	TREEP(p,j) = aux;
}

/**
 * nom : descendre
 * but : descendre un élément dans l'arbre
 * @param p : l'arbre
 * @param k : l'indice de l'élément à descendre
 * @return void
*/
void descendre(T_indirectHeap *p, int k) {
	int i, fini =  isLEAF(k,p->nbElt);
	while (!fini) {
		i = iLCHILD(k);
		if (isINTREE(iRCHILD(k), p->nbElt) && VALP(p,TREEP(p,iRCHILD(k))) > VALP(p,TREEP(p,iLCHILD(k))))
			i = iRCHILD(k);
		if (VALP(p,TREEP(p,k)) >= VALP(p,TREEP(p,i))) 
			fini = 1;
			
		else {
			swapTree(p, k, i);
			k = i;
			fini = isLEAF(k,p->nbElt); 
		}
	}
}

/**
 * nom : buildHeap
 * description : tri le tableau p->tree en utilisant la méthode du tas
 * @param p : l'arbre
 * @return void
*/
void buildHeapV2(T_indirectHeap * p){
	for (int i = p->nbElt-1; i >= 0; i--)
		descendre(p, i);

	for (int i = 0 ; i < p->nbElt; i++)
		descendre(p, i);
}

/**
 * extraireMin
 * description : extrait le minimum de l'arbre
 * @param p : l'arbre
 * @return T_elt : le minimum
*/
T_elt extraireMin(T_indirectHeap *p) {
	T_elt min = TREEP(p,0);
	swapTree(p, 0, p->nbElt-1);
	p->nbElt--;
	if (p->nbElt > 0)
		heapSortV2(p);
	return min;
}

/**
 * nom : insererMI
 * description : insère un élément dans le tree et dans le tableau de données
 * @param p : l'arbre
 * @param e : l'élément à insérer
 * @param occ : le nombre d'occurences de l'élément
 * @return void
*/
void insererMI(T_indirectHeap *p, int e, int occ){
	TREEP(p,p->nbElt) = e;
	p->nbElt++;
	VALP(p,e) = occ;
	heapSortV2(p);
}

/**
 * nom : printHuffmanTree
 * description : affiche l'arbre de codage de Huffman
 * @param p : l'arbre
 * @return void
*/
void printHuffmanTree(T_indirectHeap *p) {
	printf("Nombre d'éléments : %d\n", p->nbElt);
	for (int i = 0; i < p->nbElt; i++) {
		printf("tree %d | data(occ) %d | huff %d\n", p->tree[i], p->data[p->tree[i]], p->huffmanTree[p->tree[i]]);
	}
}

void heapSortV2(T_indirectHeap *p) {
	int taille = p->nbElt;
	buildHeapV2(p);
	while (p->nbElt > 1){
		swapTree(p, 0, p->nbElt-1);
		p->nbElt--;
		if (p->nbElt > 0)
			descendre(p, 0);
	}

	p->nbElt = taille;
}

void codageHuffman(T_indirectHeap *p) {
	//pour chaque charactère, on va créer un arbre de codage de Huffman
	//on va parcourir l'arbre de codage de Huffman pour trouver le chemin de la racine à la feuille et on va stocker ce chemin dans le tableau codage
	//dans le tableau bits, on va stocker le codage en forme de tableau de char
	int j=0,h=0,fini2=0;
	for (int i = 0; i < ILASTCAR; i++) {
		if (p->huffmanTree[i] != -256){
			h=i;
			j=1;
			fini2=0;
			while (!fini2)
			{
				if (p->huffmanTree[h]==-256){
					fini2=1;
				}
				else{
					if (p->huffmanTree[h]<0)
						p->codage[i][IMAXCODAGE-j] = '0';
					else
						p->codage[i][IMAXCODAGE-j] = '1';
					j++;
					h = abs(p->huffmanTree[h]);
				}
			}
			p->codage[i][IMAXCODAGE-j] = '\0';
		}
	}
}


void printCodage(T_indirectHeap *p) {
	printf("car : occ | long | bits\n");
	printf("----+-----+------+------\n");

	for (int i = 0; i < ILASTCAR; i++) {
		if (p->huffmanTree[i] != -256){
			printf(" %c  |   %d |      |  ", i, p->data[i]);
			for (int j = 0; j < IMAXCODAGE; j++)
				if (p->codage[i][j] != '\0')
					printf("%c", p->codage[i][j]);
			printf("\n");
		}
	}
}