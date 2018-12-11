#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <stdbool.h>
#include <sysexits.h>
#include <getopt.h>
#include <unistd.h>
#include <math.h>
#include <string.h>

typedef struct {
    int size;
    int realsize;
    int *data;
} gslist_t;

gslist_t *construct_gslist () {
    gslist_t *gslist = malloc (sizeof (gslist_t));
    gslist->size = 0;
    gslist->realsize = 1;
    gslist->data = malloc (sizeof (int));
    return gslist;
}

void destruct_gslist (gslist_t *gslist) {
    free (gslist->data);
    free (gslist);
}

void gslist_add (gslist_t *gslist, int element) {
    if (gslist->size >= gslist->realsize) {
        gslist->realsize *= 2;
        int *new = malloc (gslist->realsize * sizeof (int));
        for (int i = 0; i < gslist->size; i++)
            new[i] = gslist->data[i];
        free (gslist->data);
        gslist->data = new;
    }

    gslist->data[gslist->size++] = element;
}

typedef struct {
    int size;
    int edgecount;
    gslist_t **neighbours;
} graph_t;

graph_t *construct_graph_from_file (FILE *file) {
    char *line = NULL;
    size_t len = 0;
    ssize_t nread;

    graph_t *graph = NULL;
    int initialized = false;

    while ((nread = getline(&line, &len, file)) != -1) {
        char *tok = strtok(line, " ");

        if (strcmp (tok, "p") == 0) {
            if (initialized) {
                fprintf(stderr, "Already initialized\n");
                exit(EX_DATAERR);
            }
            initialized = true;
            strtok(NULL, " ");
            int size = atoi(strtok(NULL, " "));

            graph = malloc(sizeof(graph_t));
            graph->size = size;
            graph->neighbours = malloc(size * sizeof(gslist_t *));
            graph->edgecount = 0;
            for (int i = 0; i < size; i++)
                graph->neighbours[i] = construct_gslist();
        }

        if (strcmp (tok, "e") == 0) {
            if (!initialized) {
                fprintf(stderr, "Not yet initialized\n");
                exit(EX_DATAERR);
            }

            graph->edgecount++;

            int first = atoi(strtok(NULL, " ")) - 1;
            int second = atoi(strtok(NULL, " ")) - 1;

            gslist_add(graph->neighbours[first], second);
            gslist_add(graph->neighbours[second], first);
        }
    }

    free(line);

    if (!initialized) {
        fprintf(stderr, "Not yet initialized\n");
        exit(EX_DATAERR);
    }
    return graph;
}

void destruct_graph (graph_t *graph) {
    for (int i = 0; i < graph->size; i++)
        destruct_gslist(graph->neighbours[i]);
    free(graph->neighbours);
    free(graph);
}

typedef struct {
    int length;
    int *string;
} bitstring_t;

bitstring_t *construct_bitstring (int length) {
    bitstring_t *bitstring = malloc(sizeof(bitstring_t));
    bitstring->length = length;
    bitstring->string = malloc(length * sizeof(int));

    return bitstring;
}

bitstring_t *construct_random_bitstring (int length, int colors) {
    bitstring_t *bitstring = construct_bitstring(length);
    bitstring->length = length;
    for (int i = 0; i < length; i++) {
        bitstring->string[i] = rand() % colors;
    }

    return bitstring;
}

void destruct_bitstring (bitstring_t *bitstring) {
    free (bitstring->string);
    free (bitstring);
}

bitstring_t *copy_bitstring(bitstring_t *bitstring) {
    bitstring_t *copy = malloc(sizeof(bitstring_t));
    copy->string = malloc(bitstring->length * sizeof(int));

    for (int i = 0; i < bitstring->length; i++)
        copy->string[i] = bitstring->string[i];

    return copy;
}

int bitstring_equal(bitstring_t *first_bitstring, bitstring_t *second_bitstring) {
    // TODO: optimize this
    for (int i = 0; i < first_bitstring->length; i++)
        for (int j = 0; j < first_bitstring->length; j++) {
            if (first_bitstring->string[i] == first_bitstring->string[j]) {
                if (second_bitstring->string[i] != second_bitstring->string[j])
                    return false;
            } else {
                if (second_bitstring->string[i] == second_bitstring->string[j])
                    return false;
            }
        }

    return true;
}

typedef struct {
    int size;
    int gene_length;
    int colors;
    bitstring_t ** pool;
} genepool_t;

genepool_t *construct_genepool(int poolsize, int colors, int stringlength) {
    genepool_t *genepool = malloc(sizeof(genepool_t));
    genepool->size = poolsize;
    genepool->gene_length = stringlength;
    genepool->colors = colors;
    genepool->pool = malloc(poolsize * sizeof(bitstring_t *));

    for (int i = 0; i < poolsize; i++) {
        genepool->pool[i] = construct_random_bitstring (stringlength, colors);
    }

    return genepool;
}

genepool_t *allocate_genepool(int poolsize, int colors, int stringlength) {
    genepool_t *genepool = malloc(sizeof(genepool_t));
    genepool->gene_length = stringlength;
    genepool->size = poolsize;
    genepool->colors = colors;
    genepool->pool = malloc(poolsize * sizeof(bitstring_t *));

    for (int i = 0; i < poolsize; i++) {
        genepool->pool[i] = construct_bitstring (genepool->gene_length);
    }

    return genepool;
}

void reset_genepool(genepool_t *genepool) {
    /* Warning: This changes the bitstrings instead of replacing them.
       This saves time but is dangerous */
    for (int i = 0; i < genepool->size; i++)
        for (int j = 0; j < genepool->gene_length; j++)
            genepool->pool[i]->string[j] = random() % genepool->colors;
}

void destruct_genepool(genepool_t *genepool) {
    for (int i = 0; i < genepool->size; i++) {
        destruct_bitstring(genepool->pool[i]);
    }
    free(genepool->pool);
    free(genepool);
}

// Randomize order of all the genes
void shuffle(genepool_t *genepool) {
    // Knuth shuffle
    for (int i = 0; i < genepool->size - 1; i++) {
        // j is random between i and genepool->size-1 inclusive
        int j = i + (random() % (genepool-> size - i));
        // swap
        bitstring_t *old = genepool->pool[j];
        genepool->pool[j] = genepool->pool[i];
        genepool->pool[i] = old;
    }
}

void vertex_descent (graph_t *graph, int colors, bitstring_t *gene) {
    int improved = true;
    int *color_table = malloc(colors * sizeof(int));
    while (improved) {
        improved = false;
        for (int i = 0; i < graph->size; i++) {
            for (int j = 0; j < colors; j++)
                color_table[j] = 0;

            for (int j = 0; j < graph->neighbours[i]->size; j++)
                color_table[gene->string[graph->neighbours[i]->data[j]]]++;

            int best_value = graph->size;
            int best = 0;

            for (int j = 0; j < colors; j++)
                if (color_table[j] < best_value) {
                    best_value = color_table[j];
                    best = j;
                }

            int old_value = color_table[gene->string[i]];
            if (best_value < old_value) {
                gene->string[i] = best;
                improved = true;
            } else if (best_value == old_value && random() % 2 == 0) {
                gene->string[i] = best;
            }
        }
    }
    free(color_table);
}

int graph_coloring_fitness (graph_t *graph, bitstring_t *gene) {
    int fitness = graph->edgecount;
    for (int i = 0; i < graph->size; i++)
        for (int j = 0; j < graph->neighbours[i]->size; j++) {
            int neighbour = graph->neighbours[i]->data[j];
            if (neighbour > i && gene->string[i] == gene->string[neighbour])
                fitness--;
        }
    return fitness;
}

bitstring_t *greedy_partitioning_crossover (int colors,
                                            bitstring_t *first_parent,
                                            bitstring_t *second_parent) {
    bitstring_t *child = construct_bitstring(first_parent->length);

    // A class for each color, containing the vertices that have that color
    gslist_t **fst_classes = malloc(colors * sizeof(gslist_t *));
    gslist_t **snd_classes = malloc(colors * sizeof(gslist_t *));

    int *fst_sizes = calloc(colors, sizeof(int));
    int *snd_sizes = calloc(colors, sizeof(int));

    for (int j = 0; j < colors; j++) {
        fst_classes[j] = construct_gslist();
        snd_classes[j] = construct_gslist();
    }
    for (int j = 0; j < first_parent->length; j++) {
        int color = first_parent->string[j];
        gslist_add(fst_classes[color],j);
        fst_sizes[color]++;
        color = second_parent->string[j];
        gslist_add(snd_classes[color],j);
        snd_sizes[color]++;
    }

    int running = true;
    int side = 1;
    int color = 0;
    gslist_t **classes1, **classes2;
    int *sizes1, *sizes2;
    bitstring_t *parent2;

    while (running) {
        side = (side + 1) % 2;

        if (side == 0) {
            classes1 = fst_classes;
            classes2 = snd_classes;
            sizes1 = fst_sizes;
            sizes2 = snd_sizes;
            parent2 = second_parent;
        } else {
            classes1 = snd_classes;
            classes2 = fst_classes;
            sizes1 = snd_sizes;
            sizes2 = fst_sizes;
            parent2 = first_parent;
        }

        // Find biggest class
        int biggest = 0;
        int biggest_size = -1;
        for (int j = 0; j < colors; j++)
            if (sizes1[j] > biggest_size) {
                biggest = j;
                biggest_size = sizes1[j];
            }

        // Put biggest class into child while removing from class record
        for (int j = 0; j < classes1[biggest]->size; j++) {
            int vertex = classes1[biggest]->data[j];
            if (vertex != -1) {
                child->string[vertex] = color;

                classes1[biggest]->data[j] = -1;
                int vertexcolor2 = parent2->string[vertex];

                gslist_t *other_class = classes2[vertexcolor2];
                for (int k = 0; k < other_class->size; k++)
                    if (other_class->data[k] == vertex)
                        other_class->data[k] = -1;
                sizes2[vertexcolor2]--;
            }
        }
        sizes1[biggest] = 0;

        color++;

        running = false;
        if (color < colors)
            for (int j = 0; j < colors; j++)
                if (sizes2[j] > 0)
                    running = true;
    }

    // Leftover vertices (if color != colors we were already done)
    if (color == colors) {
        for (int i = 0; i < colors; i++)
            if (sizes1[i] > 0)
                for (int j = 0; j < classes1[i]->size; j++)
                    if (classes1[i]->data[j] != -1)
                        child->string[classes1[i]->data[j]] = random() % colors;
    }

    for (int i = 0; i < colors; i++) {
        destruct_gslist(fst_classes[i]);
        destruct_gslist(snd_classes[i]);
    }
    free(fst_classes);
    free(fst_sizes);
    free(snd_classes);
    free(snd_sizes);

    return child;
}

void crossover(genepool_t *genepool,
               graph_t *graph,
               int colors) {

    shuffle(genepool);
    for (int i = 0; i < genepool->size / 2; i++) {
        // Set parents
        bitstring_t *first_parent = genepool->pool[i*2];
        bitstring_t *second_parent = genepool->pool[i*2+1];

        if(!bitstring_equal(first_parent, second_parent)) {
            // Construct children
            bitstring_t *first_child = greedy_partitioning_crossover(colors, first_parent,second_parent);
            bitstring_t *second_child = greedy_partitioning_crossover(colors, second_parent, first_parent);

            vertex_descent(graph, colors, first_child);
            vertex_descent(graph, colors, second_child);

            int fp_fitness = graph_coloring_fitness(graph, first_parent);
            int sp_fitness = graph_coloring_fitness(graph, second_parent);
            int fc_fitness = graph_coloring_fitness(graph, first_child);
            int sc_fitness = graph_coloring_fitness(graph, second_child);

            int replace_index = 0;
            bitstring_t **replacers = malloc(2*sizeof(bitstring_t *));

            int destroy_fp = true;
            int destroy_sp = true;
            int destroy_fc = true;
            int destroy_sc = true;

            int first_chosen = 0;
            int second_chosen = 0;

            /* At the start there are 2 places available. The first child can be
               chosen if it is not strictly smaller than 2 others (equivalent to geq
               to 2 others).

               Then there are two options; either there are still two places
               available, because the first child wasn't chosen, or there is only
               one left because the first child was chosen.
               If there is only one spot left, the second child can be chosen if it
               is not strictly smaller than each parent parents.
               If there are two spots left, the second child can be chosen if it is
               not strictly smaller than at least one parent.

               Then there are three options; either there are still two places
               available, or there is one, or there are none.
               If there are still two places then the first parent is chosen.
               If there is one place left, the first parent is chosen if it is not
               strictly smaller than the second parent.
               If there are no places left, then the first parent will not be chosen.

               Then there are two options; either there is one place available, or
               there are none.
               If one place is available, the second parent is chosen.
               If no place is available, the second parent is not chosen.
            */

            if ((fc_fitness >= fp_fitness && fc_fitness >= sp_fitness) ||
                (fc_fitness >= fp_fitness && fc_fitness >= sc_fitness) ||
                (fc_fitness >= sp_fitness && fc_fitness >= sc_fitness)) {
                replacers[replace_index++] = first_child;
                destroy_fc = false;
                if (bitstring_equal(first_child,first_parent))
                    first_chosen++;
                if (bitstring_equal(first_child,second_parent))
                    second_chosen++;
            }

            if ((replace_index >= 1 &&
                 (sc_fitness >= fp_fitness && sc_fitness >= sp_fitness)) ||
                (replace_index < 1 &&
                 (sc_fitness >= fp_fitness || sc_fitness >= sp_fitness))) {
                replacers[replace_index++] = second_child;
                destroy_sc = false;
                if (bitstring_equal(second_child,first_parent))
                    first_chosen++;
                if (bitstring_equal(second_child,second_parent))
                    second_chosen++;
            }
            if (replace_index == 0 ||
                (replace_index == 1 && fp_fitness >= sp_fitness)) {
                replacers[replace_index++] = first_parent;
                destroy_fp = false;
                first_chosen++;
            }
            if (replace_index == 1) {
                replacers[replace_index++] = second_parent;
                destroy_sp = false;
                second_chosen++;
            }

            genepool->pool[i*2] = replacers[0];
            genepool->pool[i*2+1] = replacers[1];

            if (destroy_fp)
                destruct_bitstring(first_parent);
            if (destroy_sp)
                destruct_bitstring(second_parent);
            if (destroy_fc)
                destruct_bitstring(first_child);
            if (destroy_sc)
                destruct_bitstring(second_child);

            free(replacers);
        }
    }
}

float get_correlation(float pool_size, float *fitnesses_old, float *fitnesses_new) {
    float old_mean = 0, new_mean = 0;
    for (int i = 0; i < pool_size; i++) {
        old_mean += fitnesses_old[i] /  pool_size;
        new_mean += fitnesses_new[i] /  pool_size;
    }

    float cov = 0;
    float old_sd = 0;
    float new_sd = 0;

    for (int i = 0; i < pool_size; i++) {
        old_sd += pow(fitnesses_old[i] - old_mean,2) /  pool_size;
        new_sd += pow(fitnesses_new[i] - new_mean,2) /  pool_size ;
        cov += (fitnesses_old[i] - old_mean) * (fitnesses_new[i] - new_mean) / pool_size;
    }

    old_sd = pow(old_sd, 0.5);
    new_sd = pow(new_sd, 0.5);

    return cov / (old_sd * new_sd);
}

int main(int argc, char **argv) {
    srand(time(NULL) * getpid());

    char *path = NULL;
    int colors = -1;

    int get_correlations = false;

    int c;
    while ((c = getopt(argc, argv, "f:c:xh")) != -1) {
        switch (c) {
        case 'f':
            path = optarg;
            break;
        case 'c':
            colors = atoi(optarg);
            break;
        case 'x':
            get_correlations = true;
            break;
        case 'h':
            puts("WORK IN PROGRESS");
            exit(EX_OK);
            break;
        default:
            fprintf(stderr, "No handler for argument\n");
            exit(EX_USAGE);
        }
    }

    if (path == NULL) {
        fprintf(stderr, "Path to input not set\n");
        exit(EX_NOINPUT);
    }
    if (colors < 1) {
        fprintf(stderr, "Colors not valid or not set\n");
        exit(EX_NOINPUT);
    }

    FILE *file = fopen(path,"r");

    if (file == NULL) {
        fprintf(stderr, "Could not open input file \n");
        exit (EX_NOINPUT);
    }

    graph_t *graph = construct_graph_from_file(file);
    fclose(file);

    if (get_correlations) {
        int pool_size = 10000;
        genepool_t *genepool = construct_genepool(pool_size, colors, graph->size);

        float *fitnesses_old = malloc(pool_size * sizeof(int) / 2);
        float *fitnesses_new = malloc(pool_size * sizeof(int) / 2);

        // TODO is this right?
        for (int i = 0; i < pool_size / 2; i++) {
            fitnesses_old[i] =
                (float) graph_coloring_fitness(graph, genepool->pool[i*2]) *
                (float) graph_coloring_fitness(graph, genepool->pool[i*2+1]);

            bitstring_t *fst_gene = greedy_partitioning_crossover(colors,
                                                                  genepool->pool[i*2],
                                                                  genepool->pool[i*2+1]);
            bitstring_t *snd_gene = greedy_partitioning_crossover(colors,
                                                                  genepool->pool[i*2+1],
                                                                  genepool->pool[i*2]);

            fitnesses_new[i] =
                (float) graph_coloring_fitness(graph, fst_gene) *
                (float) graph_coloring_fitness(graph, snd_gene);

            destruct_bitstring(fst_gene);
            destruct_bitstring(snd_gene);
        }
        
        printf("GPX correlation: %f\n", get_correlation(pool_size/2,fitnesses_old,fitnesses_new));

        for (int i = 0; i < pool_size / 2; i++) {
            vertex_descent(graph,colors,genepool->pool[i*2]);
            vertex_descent(graph,colors,genepool->pool[i*2+1]);

            fitnesses_old[i] =
                (float) graph_coloring_fitness(graph, genepool->pool[i*2]) *
                (float) graph_coloring_fitness(graph, genepool->pool[i*2+1]);

            bitstring_t *fst_gene = greedy_partitioning_crossover(colors,
                                                                  genepool->pool[i*2],
                                                                  genepool->pool[i*2+1]);
            bitstring_t *snd_gene = greedy_partitioning_crossover(colors,
                                                                  genepool->pool[i*2+1],
                                                                  genepool->pool[i*2]);
            vertex_descent(graph,colors,fst_gene);
            vertex_descent(graph,colors,snd_gene);

            fitnesses_new[i] =
                (float) graph_coloring_fitness(graph, fst_gene) *
                (float) graph_coloring_fitness(graph, snd_gene);

            destruct_bitstring(fst_gene);
            destruct_bitstring(snd_gene);
        }


        printf("GPX+VDLS correlation: %f\n", get_correlation(pool_size/2,fitnesses_old,fitnesses_new));

        free(fitnesses_old);
        free(fitnesses_new);

        destruct_genepool(genepool);
    } else {
        int pool_size = 100;
        int found = true;
        while(found) {
            genepool_t *genepool = construct_genepool(pool_size, colors, graph->size);
            int total = 0;
            int old_total;
            int rounds_no_improvement = 0;
            int done = false;

            for (int i = 0; i < pool_size; i++)
                vertex_descent(graph,colors,genepool->pool[i]);

            while (rounds_no_improvement < 100 && !done) {
                crossover(genepool, graph, colors);

                old_total = total;
                total = 0;
                for (int i = 0; i < genepool->size; i++) {
                    int fitness = graph_coloring_fitness (graph, genepool->pool[i]);
                    total += fitness;
                    if (fitness == graph->edgecount)
                        done = true;
                }

                if(total > old_total) {
                    rounds_no_improvement = 0;
                } else {
                    rounds_no_improvement++;
                }
            }

            int max_score = graph->edgecount;
            found = false;
            for (int i = 0; i < genepool->size; i++)
                if (graph_coloring_fitness(graph,genepool->pool[i]) == max_score) {
                    found = true;
                    break;
                }

            if (found)
                printf("Found optimal solution at %i colors!\n", colors);
            else
                printf("Did not find solution at %i colors, stopping!\n", colors);
            colors--;

            destruct_genepool(genepool);
        }
    }
    destruct_graph(graph);
}
