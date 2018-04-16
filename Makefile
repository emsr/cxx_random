
all: \
  test_uniform_inside_sphere_distribution \
  test_uniform_inside_triangle_distribution \
  test_uniform_inside_tetrahedron_distribution \
  test_von_mises_fisher_distribution \
  test_dirichlet_distribution

test_uniform_inside_sphere_distribution: test_uniform_inside_sphere_distribution.cpp uniform_inside_sphere_distribution.h uniform_inside_sphere_distribution.tcc
	$(HOME)/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -o test_uniform_inside_sphere_distribution test_uniform_inside_sphere_distribution.cpp

test_uniform_inside_triangle_distribution: test_uniform_inside_triangle_distribution.cpp
	$(HOME)/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -o test_uniform_inside_triangle_distribution test_uniform_inside_triangle_distribution.cpp

test_uniform_inside_tetrahedron_distribution: test_uniform_inside_tetrahedron_distribution.cpp
	$(HOME)/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -o test_uniform_inside_tetrahedron_distribution test_uniform_inside_tetrahedron_distribution.cpp

test_von_mises_fisher_distribution: test_von_mises_fisher_distribution.cpp von_mises_fisher_distribution.h von_mises_fisher_distribution.tcc
	$(HOME)/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -o test_von_mises_fisher_distribution test_von_mises_fisher_distribution.cpp

test_dirichlet_distribution: test_dirichlet_distribution.cpp dirichlet_distribution.h dirichlet_distribution.tcc
	$(HOME)/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -o test_dirichlet_distribution test_dirichlet_distribution.cpp

test:
	./test_uniform_inside_sphere_distribution > test_uniform_inside_sphere_distribution.txt
	./test_uniform_inside_triangle_distribution > test_uniform_inside_triangle_distribution.txt
	./test_uniform_inside_tetrahedron_distribution > test_uniform_inside_tetrahedron_distribution.txt
