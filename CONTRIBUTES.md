# Contributes

Peroxide follows [Gitflow workflow](https://www.atlassian.com/git/tutorials/comparing-workflows/gitflow-workflow).
Don't need to much precise, but it's good to obey some simple rules.

1. Do not work on `master` branch directly.
2. Use `features/<features_to_improve>` branch mainly. For example, if you want to improve `DataFrame` then create `features/dataframe` on latest `dev` branch & create pull request to `dev` branch.
3. Test before push. `cargo test --features <required_features>` is wonderful.

Thanks for all contributions!
