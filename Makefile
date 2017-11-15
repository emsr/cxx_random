
all: \
  test_uniform_inside_triangle_distribution \
  test_uniform_inside_tetrahedron_distribution

test_uniform_inside_triangle_distribution: test_uniform_inside_triangle_distribution.cpp
	$(HOME)/bin/bin/g++ -std=gnu++11 -g -Wall -Wextra -o test_uniform_inside_triangle_distribution test_uniform_inside_triangle_distribution.cpp

test_uniform_inside_tetrahedron_distribution: test_uniform_inside_tetrahedron_distribution.cpp
	$(HOME)/bin/bin/g++ -std=gnu++11 -g -Wall -Wextra -o test_uniform_inside_tetrahedron_distribution test_uniform_inside_tetrahedron_distribution.cpp
