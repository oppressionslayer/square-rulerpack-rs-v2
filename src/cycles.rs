//! Live cycle counting over large primes: union-find on the graph whose vertices are large primes
//! (with 1 standing for "no prime") and whose edges are partial relations. The number of
//! independent cycles, which is the number of full relations the partials will yield, is
//! `edges + components - vertices`, kept up to date per relation so the sieve knows when to stop.
//! The same count is what YAFU and msieve print as "N from M partial".
use std::collections::HashMap;

#[derive(Default, Debug)]
pub struct CycleGraph {
    index: HashMap<u64, u32>,
    parent: Vec<u32>,
    rank: Vec<u8>,
    pub edges: u64,
    pub components: u64,
}

impl CycleGraph {
    pub fn new() -> Self { Self::default() }

    fn vertex(&mut self, p: u64) -> u32 {
        if let Some(&i) = self.index.get(&p) { return i; }
        let i = self.parent.len() as u32;
        self.index.insert(p, i);
        self.parent.push(i);
        self.rank.push(0);
        self.components += 1;
        i
    }

    fn find(&mut self, mut x: u32) -> u32 {
        let mut root = x;
        while self.parent[root as usize] != root { root = self.parent[root as usize]; }
        while self.parent[x as usize] != root { let next = self.parent[x as usize]; self.parent[x as usize] = root; x = next; }
        root
    }

    /// Add a partial relation with large primes (a, b); use 1 for a missing second prime. A full
    /// relation (1, 1) is not an edge and must not be added. Returns true if the edge closed a
    /// cycle (both primes were already connected), which is when one more full relation becomes
    /// available.
    pub fn add(&mut self, a: u64, b: u64) -> bool {
        debug_assert!(!(a == 1 && b == 1));
        let (a, b) = if a == b { (a, 1) } else { (a, b) }; // a square large prime pairs with 1
        let ia = self.vertex(a);
        let ib = self.vertex(b);
        self.edges += 1;
        let ra = self.find(ia);
        let rb = self.find(ib);
        if ra == rb { return true; }
        self.components -= 1;
        // union by rank
        match self.rank[ra as usize].cmp(&self.rank[rb as usize]) {
            std::cmp::Ordering::Less => self.parent[ra as usize] = rb,
            std::cmp::Ordering::Greater => self.parent[rb as usize] = ra,
            std::cmp::Ordering::Equal => { self.parent[rb as usize] = ra; self.rank[ra as usize] += 1; }
        }
        false
    }

    pub fn vertices(&self) -> u64 { self.parent.len() as u64 }

    /// Independent cycles so far = full relations the partials can be combined into.
    pub fn cycles(&self) -> u64 { self.edges + self.components - self.vertices() }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn single_large_prime_pairs() {
        let mut g = CycleGraph::new();
        assert!(!g.add(101, 1)); // first partial with prime 101: no cycle
        assert_eq!(g.cycles(), 0);
        assert!(g.add(101, 1)); // second with the same prime: one cycle
        assert_eq!(g.cycles(), 1);
        assert!(g.add(101, 1)); // third: another
        assert_eq!(g.cycles(), 2);
        assert!(!g.add(103, 1));
        assert_eq!(g.cycles(), 2);
    }
    #[test]
    fn double_large_prime_chain() {
        // (1,a) (a,b) (b,1): three partials, one cycle through the chain
        let mut g = CycleGraph::new();
        assert!(!g.add(1, 7)); assert!(!g.add(7, 11)); assert!(g.add(11, 1));
        assert_eq!(g.cycles(), 1);
        // a triangle a-b, b-c, c-a with no 1: one cycle
        assert!(!g.add(17, 19)); assert!(!g.add(19, 23)); assert!(g.add(23, 17));
        assert_eq!(g.cycles(), 2);
        // an edge between the two components: no new cycle
        assert!(!g.add(7, 17));
        assert_eq!(g.cycles(), 2);
        assert_eq!(g.vertices(), 6);
    }
    #[test]
    fn count_formula_matches_random_graph() {
        let mut g = CycleGraph::new();
        let mut x = 0x2545_F491_4F6C_DD1Du64;
        let mut closed = 0;
        for _ in 0..5000 {
            x ^= x << 13; x ^= x >> 7; x ^= x << 17;
            let a = 1 + x % 300; x ^= x << 13; x ^= x >> 7; x ^= x << 17;
            let b = 1 + x % 300;
            if a == 1 && b == 1 { continue; }
            if g.add(a, b) { closed += 1; }
        }
        assert_eq!(closed, g.cycles());
    }
}
