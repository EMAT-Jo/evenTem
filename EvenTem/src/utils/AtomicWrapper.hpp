// courtesy of https://stackoverflow.com/questions/13193484/how-to-declare-a-vector-of-atomic-in-c

#ifndef ATOMICWRAPPER_H
#define ATOMICWRAPPER_H

#include <atomic>
#include <vector>


template <typename T>
struct atomwrapper
{
  std::atomic<T> _a;

  atomwrapper()
    :_a()
  {}

  atomwrapper(const std::atomic<T> &a)
    :_a(a.load())
  {}

  atomwrapper(const atomwrapper &other)
    :_a(other._a.load())
  {}

  atomwrapper &operator=(const atomwrapper &other)
  {
    _a.store(other._a.load());
  }

  // Overload for ++ operator
  atomwrapper &operator++()
  {
    _a.fetch_add(1);
    return *this;
  }

  atomwrapper operator++(int)
  {
    atomwrapper temp = *this;
    ++(*this);
    return temp;
  }

  T load() const
  {
    return _a.load();
  }
};

#endif // ATOMICWRAPPER_H