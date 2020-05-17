#ifndef INCLUDED_SDSL_K2_TREE_ITERATOR
#define INCLUDED_SDSL_K2_TREE_ITERATOR

#include <tuple>
#include <deque>
#include <map>
#include <iterator>
#include <iostream>

#include <sdsl/k2_tree.hpp>

using namespace std;

namespace sdsl
{

    template <uint8_t k,
              typename t_bv = bit_vector,
              typename t_rank = typename t_bv::rank_1_type>
    class edge_iterator
    {
        using idx_type = k2_tree_ns::idx_type;
        using size_type = k2_tree_ns::size_type;
        using edge = std::tuple<idx_type, idx_type>;

    public:
        using value_type = edge;
        using pointer = edge *;
        using reference = edge &;
        using iterator_category = std::forward_iterator_tag;

        edge_iterator() {}
        edge_iterator(const t_bv &k_t, const t_bv &k_l, const t_rank &k_t_rank, const uint16_t &k_height) : k_t(&k_t), k_l(&k_l), k_t_rank(&k_t_rank), k_height(k_height)
        {
            _initialize();
        }

        value_type operator*()
        {
            return *_ptr;
        }

        edge_iterator<k, t_bv, t_rank> &operator=(const edge_iterator<k, t_bv, t_rank> &other)
        {
            if (this != &other)
            {
                this->k_t = other.k_t;
                this->k_l = other.k_l;
                this->k_t_rank = other.k_t_rank;
                this->k_height = other.k_height;

                this->_ptr = other._ptr;
                this->size = other.size;
                
                this->_state_queue = other._state_queue;
            }
            return *this;
        }

        bool operator==(const edge_iterator<k, t_bv, t_rank> &rhs) const
        {
            return equal_edge(*(rhs._ptr), *(this->_ptr));
        }
        bool operator!=(const edge_iterator<k, t_bv, t_rank> &rhs) const
        {
            return !(*this == rhs);
        }

        friend void swap(edge_iterator<k, t_bv, t_rank> &rhs, edge_iterator<k, t_bv, t_rank> &lhs)
        {
            if (lhs != rhs)
            {
                std::swap(lhs._ptr, rhs._ptr);
                std::swap(lhs.k_t, rhs.k_t);
                std::swap(lhs.k_l, rhs.k_l);
                std::swap(lhs.k_t_rank, rhs.k_t_rank);
                std::swap(lhs.size, rhs.size);

                std::swap(lhs._state_queue, rhs._state_queue);
            }
        }

        friend ostream &operator<<(ostream &os, const edge_iterator<k, t_bv, t_rank> &edg)
        {
            os << " ==== ktree Edge Iterator ==== " << endl;
            edge e = *edg._ptr;
            if (get<0>(e) == edg.size && get<1>(e) == edg.size)
            {
                os << " END NODE" << endl;
                //     os << " ============================= " << endl;
                //     return os;
            }

            os << " ptr (" << get<0>(e) << ", " << get<1>(e) << ")" << endl;
            os << " CONTAINER " << endl;
            os << "     k_t ";
            if (edg.k_t != NULL)
            {
                for (uint i = 0; i < edg.k_t->size(); i++)
                    os << (*edg.k_t)[i];
                os << endl;
            }
            else
            {
                os << "null" << endl;
            }
            os << "     k_l ";
            if (edg.k_l != NULL)
            {
                for (uint i = 0; i < edg.k_l->size(); i++)
                    os << (*edg.k_l)[i];
                os << endl;
            }
            else
            {
                os << "null" << endl;
            }
            os << "     k_height " << edg.k_height << endl;
            os << " ============================= " << endl;
            return os;
        }

        edge_iterator<k, t_bv, t_rank> end()
        {
            edge_iterator<k, t_bv, t_rank> it = *this;
            it._ptr = new edge(size, size); //end node
            return it;
        }

        edge_iterator<k, t_bv, t_rank> &operator++(int)
        {
            edge_iterator<k, t_bv, t_rank> *tmp;
            tmp = new edge_iterator<k, t_bv, t_rank>(*(this->k_t), *(this->k_l), *(this->k_t_rank), this->k_height);
            operator++();
            return *tmp;
        }

        edge_iterator<k, t_bv, t_rank> &operator++()
        {
            while (!_state_queue.empty())
            {
                state st = _state_queue.front();
                _state_queue.pop_front();

                idx_type neigh = size;
                _find_next_recursive(st.n, st.row, st.col, st.level, neigh, st.j);
                if (neigh != size)
                {
                    _ptr = new edge(st.node, neigh);
                    return *this;
                }
            }
            if (_state_queue.empty())
                _ptr = new edge(size, size);
            return *this;
        }

    protected:
        typedef struct state
        {
            idx_type node;
            unsigned row, col, n, level, j, y;

            void set(idx_type _node, unsigned _n, unsigned _row, unsigned _col, unsigned _level,
                     unsigned _j, unsigned _y)
            {
                node = _node;
                row = _row;
                col = _col;
                n = _n;
                level = _level;
                j = _j;
                y = _y;
            }

            void set(idx_type _node, unsigned _n, unsigned _row, unsigned _col, unsigned _level)
            {
                node = _node;
                row = _row;
                col = _col;
                n = _n;
                level = _level;
            }

            bool operator==(const state &s) const
            {
                return this->node == s.node && this->n == s.n && this->row == s.row && this->col == s.col;
            }

            bool operator!=(const state &s) const
            {
                return !(*this == s);
            }
        } state;

        struct stateCmp
        {
            bool operator()(const state &lhs, const state &rhs) const
            {
                unsigned a = lhs.node+1 + lhs.n + lhs.row + lhs.col + lhs.level *2;
                unsigned b = rhs.node+1 + rhs.n + rhs.row + rhs.col + rhs.level *2;
                return a < b;
            }
        };
        // container
        const t_bv *k_t = NULL;
        const t_bv *k_l = NULL;
        const t_rank *k_t_rank = NULL;
        uint16_t k_height;
        //
        // iterator state //
        pointer _ptr;
        deque<state> _state_queue;
        size_t size;
        //

        void
        _initialize()
        {
            size = std::pow(k, k_height);
            _ptr = new edge(size, size);

            if (k_l->size() > 0)
            {
                _find_next_queue();
                operator++();           
            }
        }

        void _find_next_queue()
        {
            map<state, bool, struct stateCmp> _duplicate_st;
            idx_type node = 0;
            unsigned row = 0;
            unsigned col = 0;
            unsigned n = static_cast<size_type>(std::pow(k, k_height)) / k;
            unsigned level = k * std::floor(node / static_cast<double>(n));

            for (; node < size; row = 0, col = 0, node++, level = k * std::floor(node / static_cast<double>(n)))
                for (; row < k; row++)
                    for (col = 0; col < k; col++)
                        _find_next_recursive_queue(n / k, node % n, n * row, level + row, col, node, _duplicate_st);
        }

        void _find_next_recursive_queue(unsigned n, unsigned row, unsigned col, unsigned level, unsigned initial_j, idx_type node, map<state, bool, struct stateCmp> &_duplicate_st)
        {
            if (level >= k_t->size()) // Last level
                return;

            if ((*k_t)[level] == 1)
            {
                size_type y = (*k_t_rank)(level + 1) * k * k +
                              k * std::floor(row / static_cast<double>(n));
                for (unsigned j = initial_j; j < k; j++)
                {
                    state st;
                    st.set(node, n / k, row % n, col + n * j, y + j);
                    if (n == 1 && !_duplicate_st[st])
                    {
                        _duplicate_st[st] = true;
                        _state_queue.push_back(st);
                    }
                    else
                        _find_next_recursive_queue(n / k, row % n, col + n * j, y + j, 0, node, _duplicate_st);
                }
            }
        }

        bool _find_next_recursive(unsigned n, unsigned row, unsigned col, unsigned level, idx_type &neigh, unsigned initial_j)
        {
            if (level >= k_t->size()) // Last level
            {
                if ((*k_l)[level - k_t->size()] == 1)
                {
                    neigh = col;
                    return true;
                }
                return false;
            }

            if ((*k_t)[level] == 1)
            {
                unsigned y = (*k_t_rank)(level + 1) * k * k +
                              k * std::floor(row / static_cast<double>(n));
                for (unsigned j = initial_j; j < k; j++)
                    if (_find_next_recursive(n / k, row % n, col + n * j, y + j, neigh, 0))
                        return true;
            }
            return false;
        }

    private:
        bool equal_edge(const edge &e1, const edge &e2) const
        {
            idx_type e1_x = std::get<0>(e1);
            idx_type e1_y = std::get<1>(e1);

            idx_type e2_x = std::get<0>(e2);
            idx_type e2_y = std::get<1>(e2);

            return e1_x == e2_x && e1_y == e2_y;
        }
    };
} // namespace sdsl

#endif