import { createAndPatchRegistry } from "@mat3ra/code/dist/js/context";

import { wodeProviders } from "./providers";

export const ContextProviderRegistry = createAndPatchRegistry(wodeProviders);
