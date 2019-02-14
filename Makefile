
INC_DIR = include/ext

TEST_BIN_DIR = bin

TEST_OUT_DIR = test_output

all: $(TEST_BIN_DIR) \
  $(TEST_BIN_DIR)/test_uniform_inside_sphere_distribution \
  $(TEST_BIN_DIR)/test_uniform_inside_triangle_distribution \
  $(TEST_BIN_DIR)/test_uniform_inside_tetrahedron_distribution \
  $(TEST_BIN_DIR)/test_von_mises_fisher_distribution \
  $(TEST_BIN_DIR)/test_dirichlet_distribution

$(TEST_BIN_DIR)/test_uniform_inside_sphere_distribution: test_uniform_inside_sphere_distribution.cpp $(INC_DIR)/uniform_inside_sphere_distribution.h $(INC_DIR)/uniform_inside_sphere_distribution.tcc
	$(HOME)/bin/bin/g++ -std=gnu++17 -Iinclude -g -Wall -Wextra -o $(TEST_BIN_DIR)/test_uniform_inside_sphere_distribution test_uniform_inside_sphere_distribution.cpp

$(TEST_BIN_DIR)/test_uniform_inside_triangle_distribution: test_uniform_inside_triangle_distribution.cpp $(INC_DIR)/uniform_inside_triangle_distribution.h
	$(HOME)/bin/bin/g++ -std=gnu++17 -Iinclude -g -Wall -Wextra -o $(TEST_BIN_DIR)/test_uniform_inside_triangle_distribution test_uniform_inside_triangle_distribution.cpp

$(TEST_BIN_DIR)/test_uniform_inside_tetrahedron_distribution: test_uniform_inside_tetrahedron_distribution.cpp $(INC_DIR)/uniform_inside_tetrahedron_distribution.h
	$(HOME)/bin/bin/g++ -std=gnu++17 -Iinclude -g -Wall -Wextra -o $(TEST_BIN_DIR)/test_uniform_inside_tetrahedron_distribution test_uniform_inside_tetrahedron_distribution.cpp

$(TEST_BIN_DIR)/test_von_mises_fisher_distribution: test_von_mises_fisher_distribution.cpp $(INC_DIR)/von_mises_fisher_distribution.h $(INC_DIR)/von_mises_fisher_distribution.tcc
	$(HOME)/bin/bin/g++ -std=gnu++17 -Iinclude -g -Wall -Wextra -o $(TEST_BIN_DIR)/test_von_mises_fisher_distribution test_von_mises_fisher_distribution.cpp

$(TEST_BIN_DIR)/test_dirichlet_distribution: test_dirichlet_distribution.cpp $(INC_DIR)/dirichlet_distribution.h $(INC_DIR)/dirichlet_distribution.tcc
	$(HOME)/bin/bin/g++ -std=gnu++17 -Iinclude -g -Wall -Wextra -o $(TEST_BIN_DIR)/test_dirichlet_distribution test_dirichlet_distribution.cpp

test: $(TEST_OUT_DIR)
	./$(TEST_BIN_DIR)/test_uniform_inside_sphere_distribution > $(TEST_OUT_DIR)/test_uniform_inside_sphere_distribution.txt
	./$(TEST_BIN_DIR)/test_uniform_inside_triangle_distribution > $(TEST_OUT_DIR)/test_uniform_inside_triangle_distribution.txt
	./$(TEST_BIN_DIR)/test_uniform_inside_tetrahedron_distribution > $(TEST_OUT_DIR)/test_uniform_inside_tetrahedron_distribution.txt
	./$(TEST_BIN_DIR)/test_von_mises_fisher_distribution > $(TEST_OUT_DIR)/test_von_mises_fisher_distribution.txt
	./$(TEST_BIN_DIR)/test_dirichlet_distribution > $(TEST_OUT_DIR)/test_dirichlet_distribution.txt

clean:
	rm -rf $(TEST_BIN_DIR)/*

$(TEST_BIN_DIR):
	if test ! -d $(TEST_BIN_DIR); then \
	  mkdir $(TEST_BIN_DIR); \
	fi

$(TEST_OUT_DIR):
	if test ! -d $(TEST_OUT_DIR); then \
	  mkdir $(TEST_OUT_DIR); \
	fi
