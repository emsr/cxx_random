
all: \
  test_uniform_inside_sphere_distribution \
  test_uniform_inside_triangle_distribution \
  test_uniform_inside_tetrahedron_distribution

test_uniform_inside_sphere_distribution: test_uniform_inside_sphere_distribution.cpp
	$(HOME)/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -o test_uniform_inside_sphere_distribution test_uniform_inside_sphere_distribution.cpp

test_uniform_inside_triangle_distribution: test_uniform_inside_triangle_distribution.cpp
	$(HOME)/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -o test_uniform_inside_triangle_distribution test_uniform_inside_triangle_distribution.cpp

test_uniform_inside_tetrahedron_distribution: test_uniform_inside_tetrahedron_distribution.cpp
	$(HOME)/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -o test_uniform_inside_tetrahedron_distribution test_uniform_inside_tetrahedron_distribution.cpp

test:
	./test_uniform_inside_sphere_distribution > test_uniform_inside_sphere_distribution.txt
	./test_uniform_inside_triangle_distribution > test_uniform_inside_triangle_distribution.txt
	./test_uniform_inside_tetrahedron_distribution > test_uniform_inside_tetrahedron_distribution.txt
