enum SumcheckSystem {
    Single,
    Multi,
    Cuda,
}

impl SumcheckSystem {
    fn all() -> Vec<SumcheckSystem> {
        vec![
            SumcheckSystem::Single,
            SumcheckSystem::Multi,
            SumcheckSystem::Cuda,
        ]
    }

    fn bench(&self, k: usize) {}
}
