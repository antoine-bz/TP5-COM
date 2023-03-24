#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <libgen.h> // basename
#include <math.h> 


#define MAXCARS 128
#define MAXLENGTH 200000
#define ILASTCAR MAXCARS-1
#define IMAXCODAGE 16

//Macro-fonctions de debug pour graphviz
#ifndef DEBUG
#define CHECK_IF(stat, val, msg)    \
	if ( (stat) == (val) )          \
	{                               \
		perror(msg);                \
		exit(EXIT_FAILURE);       \
    }                       
#else
#define CHECK_IF(stat, val, msg)    \
	if ( (stat) == (val) )          \
	{                               \
    	perror(msg);                \
		exit(EXIT_FAILURE);       \
    }                               \
    else printf("%s ... OK\n", msg)
#endif

typedef char T_octet[IMAXCODAGE];

typedef struct {
	int longueur;
	T_octet codage;
} T_codage;

typedef struct {
	unsigned int nbElt;
	unsigned char tree[MAXCARS]; //tableau des caractères
	int data[2*MAXCARS-1]; //tableau des occurences
	int huffmanTree[2*MAXCARS -1];
	T_codage codage[MAXCARS]; //tableau des codages binaire
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
#define VALP(heap, i)		heap->data[i]

//Prototypes
static void genDotPOThuff_rec(T_indirectHeap *p, int root, FILE *fp); 
static void genDotPOTtree_rec(T_indirectHeap *p, int root, FILE *fp);
void createDotPOT(T_indirectHeap *p,const char *basename, int huffmanTree); 


T_indirectHeap *analyserDocument(char *document);
void huffman(char *document, char *output);

T_indirectHeap * newHeap();
void swapTree(T_indirectHeap *p, int i, int j);
void descendre(T_indirectHeap *p, int k);
void buildHeapV2(T_indirectHeap *p);//transformerEnMinimierV2
void heapSortV2(T_indirectHeap *p);
unsigned char extraireMin(T_indirectHeap *p);
void insererMI(T_indirectHeap *p, int e, int occ); 

void codageHuffman(T_indirectHeap *p);
void printCodage(T_indirectHeap *p);

void encodageDocument(T_indirectHeap *p, char *document, char *buffer);
int lenFile(char *document);
void ecrireDocument(char filename[FILENAME_MAX], char *document);
void compressionDocument(char filename[FILENAME_MAX], char output[FILENAME_MAX]);
void decompressionDocument(char filename[FILENAME_MAX]);


char * outputPath = "./tp5";



int main(int argc, char *argv[]) {
	if (argc == 3 || argc == 2) {
		//on verifie si le fichier source existe
		FILE *f = fopen(argv[1], "r");
		if (f == NULL) {
			printf("Le fichier source n'existe pas\n");
			return 0;
		}
		else if (argc == 3) {
			compressionDocument(argv[1], argv[2]);
		}
		else {
			decompressionDocument(argv[1]);
		}
		fclose(f);
		
	}
	// printf("out: %s\n", buffer);
	return 0;
}



/**
 * nom : genDotPOTtree_rec
 * description : génère le fichier dot correspondant à l'arbre de codage de Huffman
 * @param p : le tas
 * @param root : indice de la racine du sous-arbre à produire
 * @param fp : flux correspondant à un fichier ouvert en écriture où écrire le sous-arbre
*/
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

/**
 * nom : genDotPOTtree_rec
 * description : génère le fichier dot correspondant aux lettres et à leurs fréquences
 * @param p : le tas
 * @param root : indice de la racine du sous-arbre à produire
 * @param fp : flux correspondant à un fichier ouvert en écriture où écrire le sous-arbre
*/
static void genDotPOTtree_rec(T_indirectHeap *p, int root, FILE *fp){
	// Attention : les fonction toString utilisent un buffer alloué comme une variable statique 
	// => elles renvoient toujours la même adresse 
	// => on ne peut pas faire deux appels à toString dans le même printf()
	
	// t : tas
	// n : taille du tas
	// root : indice de la racine du sous-arbre à produire
	// fp : flux correspondant à un fichier ouvert en écriture où écrire le sous-arbre

	// ordre de récurrence  ?

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

/**
 * nom : createDotPOT
 * description : génère le fichier dot ET PNG correspondant aux arbres
 * @param p : le tas
 * @param basename : nom de base du fichier dot et du fichier png
 * @param huffmanTree : 1 si on veut générer le fichier dot correspondant à l'arbre de Huffman, 0 sinon
 * @return : void
*/
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
 * @param document : le document à coder
 * @param output : le fichier de sortie
 * @return : void
*/
void huffman(char *document, char *output) {
	int C1, C2, n, Ni;

	T_indirectHeap *Mi = analyserDocument(document);
	heapSortV2(Mi);

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
	}
	
	createDotPOT(Mi, "huffmanTree",1);
	//Creation du tableau de codage
	codageHuffman(Mi);
	//Affichage du tableau de codage
	printCodage(Mi);
	//Decodage du document
	encodageDocument(Mi, document,output);

	printf("Longueur du code binaire : %d bits \n", lenFile(document)*8);
	printf("Longueur du code de huffman : %ld bits \n", strlen(output));
	printf("Ratio de compression : %.2f%%\n", (float)strlen(output)/(lenFile(document)*8)*100);
	free(Mi);
}

/**
 * nom : analyserDocument
 * but : compter le nombre d'occurences de chaque caractère dans le document et les stocker les données dans les tableaux
 * @param document : nom du fichier à analyse
 * @return : T_indirectHeap *Mi
*/
T_indirectHeap *analyserDocument(char *document) {
	FILE *f=fopen(document,"rt");
    char s[MAXLENGTH];
	T_indirectHeap *Mi = newHeap();
	int trouve = 0, intC;
    
    while(!feof(f)){
		fgets(s, MAXLENGTH, f );
        for (int i = 0; i < strlen(s) ; i++) {			
			trouve = 0;
			//On convertit le caractère en entier
			intC = (int)s[i];

			//On vérifie si le caractère est déjà présent dans le tableau
			for (int j = 0; j < Mi->nbElt; j++) {
				if (Mi->tree[j] == intC) {
					Mi->data[intC]++;
					trouve = 1;
					break;
				}
			}

			//Si le caractère n'est pas présent dans le tableau, on l'ajoute
			if (trouve == 0) {
				Mi->tree[Mi->nbElt] = intC;
				Mi->data[intC] = 1;
				Mi->nbElt++;
			}
		}
    }
    fclose(f);
	return Mi;
}

/**
 * nom : newHeap
 * description : alloue un tas
 * @return : le tas alloué
 */
T_indirectHeap * newHeap(){
	T_indirectHeap * pAux= (T_indirectHeap *) malloc(sizeof(T_indirectHeap));
	pAux->nbElt = 0;
	for (int i = 0; i < 2*MAXCARS -1; i++) {
		pAux->huffmanTree[i] = -256;
	}
	return pAux;
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
	unsigned char aux = TREEP(p,i);
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
}

/**
 * extraireMin
 * description : extrait le minimum de l'arbre
 * @param p : l'arbre
 * @return unsigned char : le minimum
*/
unsigned char extraireMin(T_indirectHeap *p) {
	unsigned char min = TREEP(p,0);
	swapTree(p, 0, p->nbElt-1);
	p->nbElt--;
	if (p->nbElt > 0)
		heapSortV2(p);
	return min;
}

/**
 * nom : heapSort
 * description : tri le tableau p->tree en utilisant la méthode du tas
 * @param p : l'arbre
 * @return void
*/
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
 * nom : codageHuffman
 * description : créer le tableau de codage de Huffman
 * @param p : l'arbre
 * @return void
*/
void codageHuffman(T_indirectHeap *p) {
	//pour chaque charactère, on va créer un arbre de codage de Huffman
	//on va parcourir l'arbre de codage de Huffman pour trouver le chemin de la racine à la feuille et on va stocker ce chemin dans le tableau codage
	//dans le tableau bits, on va stocker le codage en forme de tableau de char
	int j=0,h=0,l=0,fini2=0;
	for (int i = 0; i < ILASTCAR; i++) {
		if (p->huffmanTree[i] != -256){
			h=i;
			j=1;
			fini2=0;
			p->codage[i].longueur = 0;
			while (!fini2)
			{
				if (p->huffmanTree[h]==-256){
					fini2=1;
				}
				else{
					if (p->huffmanTree[h]<0)
						p->codage[i].codage[IMAXCODAGE-j] = '0';
					else
						p->codage[i].codage[IMAXCODAGE-j] = '1';				
					j++;
					h = abs(p->huffmanTree[h]);
				}
			}
			p->codage[i].longueur= j-1;	
			p->codage[i].codage[IMAXCODAGE-j] = '\0';

			l=0;
			for (int k = p->codage[i].longueur; k > 0; k--){
				p->codage[i].codage[l] = p->codage[i].codage[IMAXCODAGE-k];
				p->codage[i].codage[IMAXCODAGE-k] = '\0';
				l++;
			}
				
		}
	}
}

/**
 * nom : printCodage
 * description : affiche le tableau de codage de Huffman
 * @param p : l'arbre
 * @return void
*/
void printCodage(T_indirectHeap *p) {
	printf("\ncar : occ | long | bits\n");
	printf("----+-----+------+------\n");

	for (int i = 0; i < ILASTCAR; i++) {
		if (p->huffmanTree[i] != -256){
			printf(" %c  |   %d |   %d  | %s", i, p->data[i], p->codage[i].longueur,p->codage[i].codage);
			
			printf("\n");
		}
	}
	printf("\n");
}

/**
 * nom : encodageDocument
 * description : décode un document
 * @param p : l'arbre pour le codage
 * @param document : le document à décoder
 * @param buffer : le buffer dans lequel on va stocker le document décodé
 * @return void
*/
void encodageDocument(T_indirectHeap *p, char *document,char *buffer){
	int car,j;
	FILE *f=fopen(document,"rt");
    char s[MAXLENGTH];

	j=0;
	while (!feof(f))
	{
		fgets(s,MAXLENGTH,f);
		for (int i = 0; i < strlen(s); i++){
			car = s[i];
			for (int k = 0; k < p->codage[car].longueur; k++){
				buffer[j] = p->codage[car].codage[k];
				j++;
			}
		}
	}
	fclose(f);		
}

/**
 * nom : lenFile
 * description : calcule le nombre de caractères d'un fichier
 * @param document : le document à calculer
 * @return int : la longueur du fichier
*/
int lenFile(char *document){
	int j;
	FILE *f=fopen(document,"rt");
	char s[MAXLENGTH];

	j=0;
	while (!feof(f))
	{
		fgets(s,MAXLENGTH,f);
		for (int i = 0; i < strlen(s); i++)j++;
	}
	fclose(f);	
	return j;
}

/**
 * nom : ecrireDocument
 * description : compresser dans un fichier le document
 * @param filename : chemin du fichier à produire
 * @param document : le document compressé
 * @return void
*/
void ecrireDocument(char filename[FILENAME_MAX], char *document){
	// on verifie si le fichier existe sinon on le crée
	FILE *f = fopen(filename, "wb");
	 int i;
    for (i = 0; i < strlen(document); i += 8) {  // parcours la chaîne de codage par groupes de 8 caractères (octets)
        char byte_str[9] = {0};  // chaîne de caractères pour stocker un octet binaire
        int j;
        for (j = 0; j < 8; j++) {
            byte_str[j] = document[i+j];  // copie les 8 caractères dans la chaîne de l'octet
        }
        unsigned char byte = strtol(byte_str, NULL, 2);  // conversion de la chaîne binaire en octet
        fwrite(&byte, sizeof(unsigned char), 1, f);  // écriture de l'octet dans le fichier binaire
    }
    
	fclose(f);
}

/**
 * nom : compressionDocument
 * description : compresser un fichier
 * @param filename : chemin du fichier à compresser
 * @param output : chemin du fichier à produire
 * @return void
*/
void compressionDocument(char filename[FILENAME_MAX], char output[FILENAME_MAX]){
	char * buffer = (char *) malloc(sizeof(char) * MAXLENGTH);
	huffman(filename, buffer);
	ecrireDocument(output, buffer);
}

/**
 * nom : decompressionDocument
 * description : décompresser un fichier
 * @param filename : chemin du fichier à décompresser
 * @return void
*/
void decompressionDocument(char filename[FILENAME_MAX]){
	FILE *f=fopen(filename,"rb");
    unsigned char byte;
    char binary[9]; // chaine de caractère pour stocker la représentation binaire de chaque byte
    binary[8] = '\0'; // ajouter un caractère de fin de chaine
    char* binary_str = malloc(sizeof(char) * (lenFile(filename) * 8 + 1)); // allouer de la mémoire pour stocker la chaine de caractère binaire
    binary_str[0] = '\0'; // initialiser la chaine de caractère binaire

    while (fread(&byte, sizeof(byte), 1, f) == 1) {
        for (int i = 7; i >= 0; i--) {
            binary[i] = (byte & 1) ? '1' : '0'; // convertir chaque bit en caractère binaire ('0' ou '1')
            byte >>= 1; // décaler le byte de 1 bit vers la droite pour traiter le bit suivant
        }
        strcat(binary_str, binary); // ajouter la représentation binaire du byte à la chaine de caractère binaire complète
    }

    printf("%s\n", binary_str); // afficher la chaine de caractère binaire complète

    fclose(f);
    free(binary_str);
}
